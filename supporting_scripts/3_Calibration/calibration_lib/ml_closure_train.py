# Copyright 2026, Diogo Costa
# This file is part of OpenWQ model.
#
# openWQ hybrid physics-ML — LAYER 2 (learned flux closures), TRAINING side.
"""
Differentiable emulator + trainer for a Layer-2 flux closure.

A learned closure corrects a physics flux multiplicatively and BOUNDEDLY::

    flux_corrected = Phi_physics * (1 + alpha * g(inputs))          alpha=0 => physics

To train ``g`` end-to-end against observations, the physics has to be
DIFFERENTIABLE: we unroll a small reactive system in time with the corrected
flux and back-propagate the trajectory error through the simulation (BPTT) into
``g``'s weights. openWQ itself never differentiates — training happens here in
Python; the trained weights are exported (``export_mlp_weights``) to the
``ML_CLOSURE`` config that openWQ's ``OpenWQ_ML_closure`` reads for inference.

This module is framework-free (numpy only, manual gradients) so it runs anywhere
the calibration stack runs. The gradients are finite-difference checked
(``gradient_check``) so the hand-derived BPTT is trustworthy.

The emulator here is a 1-box first-order reactive system (concentration decay),
which is enough to demonstrate + validate the pipeline; the SAME trainer shape
(rollout -> loss -> BPTT -> export) generalizes to a real module's flux by
swapping the ``step`` function for that module's (differentiable) flux.
"""

import numpy as np

try:                                   # works as a package module or standalone
    from . import ml_regionalization
except ImportError:                    # pragma: no cover
    import ml_regionalization


# ===========================================================================
# Closure network: g(x) = tanh(W2 · tanh(W1·x + b1) + b2)  -> bounded in [-1, 1]
# (Exported as two tanh layers, exactly what OpenWQ_ML evaluates.)
# ===========================================================================
class ClosureMLP:
    def __init__(self, n_in=1, n_hidden=8, seed=0):
        rng = np.random.default_rng(seed)
        s = 0.5
        self.W1 = rng.normal(0.0, s, (n_hidden, n_in))
        self.b1 = np.zeros(n_hidden)
        self.W2 = rng.normal(0.0, s, (1, n_hidden))
        self.b2 = np.zeros(1)

    # ---- forward: single input vector x -> (g scalar, cache) ----
    def forward(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        z1 = self.W1 @ x + self.b1
        h1 = np.tanh(z1)
        z2 = float(self.W2 @ h1 + self.b2)
        g = np.tanh(z2)
        return g, (x, z1, h1, z2, g)

    # ---- backward: upstream dL/dg -> (weight grads, dL/dx) ----
    def backward(self, cache, dL_dg):
        x, z1, h1, z2, g = cache
        dz2 = dL_dg * (1.0 - g * g)                 # tanh'
        dW2 = (dz2 * h1).reshape(1, -1)
        db2 = np.array([dz2])
        dh1 = self.W2.reshape(-1) * dz2
        dz1 = dh1 * (1.0 - h1 * h1)
        dW1 = np.outer(dz1, x)
        db1 = dz1
        dL_dx = self.W1.T @ dz1                      # gradient w.r.t. the input
        return {"W1": dW1, "b1": db1, "W2": dW2, "b2": db2}, dL_dx

    # ---- parameter helpers ----
    def params(self):
        return {"W1": self.W1, "b1": self.b1, "W2": self.W2, "b2": self.b2}

    def step_sgd(self, grads, lr):
        for k in ("W1", "b1", "W2", "b2"):
            getattr(self, k)[...] -= lr * grads[k]

    def export_weights(self):
        return ml_regionalization.export_mlp_weights(
            [(self.W1, self.b1, "tanh"), (self.W2, self.b2, "tanh")])

    def to_ml_closure(self, alpha, max_correction=1.0, weights_file=None):
        """Build the openWQ ML_CLOSURE config block.

        For a REAL openWQ config pass ``weights_file``: the trained network is
        written there (JSON) and the block references it via ``WEIGHTS_FILEPATH``
        — openWQ's config normalizer upper-cases keys/values and cannot hold the
        nested weights inline, but leaves ``*FILEPATH`` values untouched. Without
        ``weights_file`` the weights are returned inline (for direct /
        non-normalized use such as unit tests). ``alpha=0`` => pure physics
        (there is no on/off flag)."""
        weights = self.export_weights()
        if weights_file:
            import json as _json
            with open(weights_file, "w") as _f:
                _json.dump(weights, _f)
            return {"ALPHA": float(alpha), "MAX_CORRECTION": float(max_correction),
                    "WEIGHTS_FILEPATH": weights_file}
        return {"alpha": float(alpha), "max_correction": float(max_correction),
                "weights": weights}


# ===========================================================================
# Pluggable differentiable physics step
# ---------------------------------------------------------------------------
# A `PhysicsStep` advances ONE state variable by one timestep under the
# corrected flux and exposes the two local derivatives the BPTT needs:
#   forward(C, g) -> C_next
#   grads(C, g)   -> (dC_next/dC  [NOT through g],  dC_next/dg)
# Swap this class for a specific module's (differentiable) flux to train a
# production closure; the trainer below is otherwise flux-agnostic. The closure
# input is still the scalar state C here (1 feature) for simplicity — extend to
# a state+attribute vector by widening ClosureMLP's n_in and the step signature.
# ===========================================================================
class PhysicsStep:
    """Interface (see FirstOrderDecay for the reference implementation)."""
    def forward(self, C, g):            # pragma: no cover
        raise NotImplementedError
    def grads(self, C, g):              # pragma: no cover
        raise NotImplementedError


class FirstOrderDecay(PhysicsStep):
    """C_{t+1} = C_t * (1 - dt*k*(1 + alpha*g)), clamped at 0. First-order decay
    corrected by the closure — the reference emulator for validation."""
    def __init__(self, k, alpha, dt):
        self.k, self.alpha, self.dt = float(k), float(alpha), float(dt)

    def forward(self, C, g):
        Cn = C * (1.0 - self.dt * self.k * (1.0 + self.alpha * g))
        return Cn if Cn > 0.0 else 0.0

    def grads(self, C, g):
        d_dC = 1.0 - self.dt * self.k * (1.0 + self.alpha * g)   # direct
        d_dg = -self.dt * self.k * self.alpha * C                # through g
        return d_dC, d_dg


def rollout(mlp, C0, n_steps, step_model):
    """Forward the corrected physics with a pluggable step. Returns
    (trajectory C[t], per-step caches)."""
    C = np.zeros(n_steps + 1)
    C[0] = C0
    caches = []
    for t in range(n_steps):
        g, cache = mlp.forward([C[t]])
        caches.append((g, cache))
        C[t + 1] = step_model.forward(C[t], g)
    return C, caches


def loss_and_grads(mlp, C0, obs, step_model):
    """MSE of the corrected trajectory vs obs + BPTT gradients w.r.t. g's weights.
    `step_model` is any PhysicsStep (default emulator: FirstOrderDecay)."""
    n_steps = len(obs) - 1
    C, caches = rollout(mlp, C0, n_steps, step_model)
    N = len(obs)
    loss = float(np.mean((C - obs) ** 2))

    grads = {kk: np.zeros_like(v) for kk, v in mlp.params().items()}
    lam = 0.0  # dL/dC[t+1] carried backward
    for t in range(n_steps, -1, -1):
        lam += 2.0 * (C[t] - obs[t]) / N            # direct loss term at step t
        if t == 0:
            break
        g, cache = caches[t - 1]                     # g used to make C[t] from C[t-1]
        d_dC, d_dg = step_model.grads(C[t - 1], g)
        dL_dg = lam * d_dg
        gw, dL_dCprev_via_g = mlp.backward(cache, dL_dg)
        for kk in grads:
            grads[kk] += gw[kk]
        # the closure input is the scalar C[t-1] (1 feature) -> collapse the
        # input-gradient vector to a scalar.
        lam = lam * d_dC + float(dL_dCprev_via_g[0])
    return loss, grads, C


def gradient_check(seed=0, eps=1e-6, step_model=None):
    """Finite-difference check of the BPTT gradients (returns max rel error)."""
    mlp = ClosureMLP(n_in=1, n_hidden=5, seed=seed)
    rng = np.random.default_rng(seed + 1)
    obs = np.abs(rng.normal(5.0, 1.0, 12))
    if step_model is None:
        step_model = FirstOrderDecay(k=0.3, alpha=0.5, dt=0.1)
    C0 = 6.0
    _, grads, _ = loss_and_grads(mlp, C0, obs, step_model)
    max_rel = 0.0
    for name in ("W1", "b1", "W2", "b2"):
        P = getattr(mlp, name)
        it = np.nditer(P, flags=["multi_index"])
        while not it.finished:
            idx = it.multi_index
            orig = P[idx]
            P[idx] = orig + eps
            lp, _, _ = loss_and_grads(mlp, C0, obs, step_model)
            P[idx] = orig - eps
            lm, _, _ = loss_and_grads(mlp, C0, obs, step_model)
            P[idx] = orig
            num = (lp - lm) / (2 * eps)
            ana = grads[name][idx]
            denom = max(1e-8, abs(num) + abs(ana))
            max_rel = max(max_rel, abs(num - ana) / denom)
            it.iternext()
    return max_rel


def train(mlp, C0, obs, step_model, iters=2000, lr=0.02, verbose=False):
    """Gradient-descent training of the closure to fit the observed trajectory."""
    hist = []
    for i in range(iters):
        loss, grads, _ = loss_and_grads(mlp, C0, obs, step_model)
        mlp.step_sgd(grads, lr)
        hist.append(loss)
        if verbose and (i % max(1, iters // 10) == 0):
            print(f"  iter {i:5d}  loss {loss:.6e}")
    return hist
