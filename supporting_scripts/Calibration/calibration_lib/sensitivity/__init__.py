# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

"""
Sensitivity Analysis
====================

Morris and Sobol sensitivity analysis for parameter screening.
"""

from .morris_screening import MorrisScreening, MorrisResult
from .sobol_analysis import SobolAnalysis, SobolResult
