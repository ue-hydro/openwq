#include "models_CH/headerfile_CH.hpp"
#include "IPhreeqc.hpp"
#include <cmath>  // for isnan, isinf
#include <cstdio>  // fprintf (OPENWQ_DEBUG_VOLR)
#include <cstdlib> // getenv
#include <algorithm>
#include <cctype>
#include <map>
#include <unistd.h>  // for getcwd
#include <climits>   // for PATH_MAX

// Host dependency names are matched without regard to case: the module JSON is
// upper-cased when read ("Tsoil_K" arrives as "TSOIL_K"), the host registers
// its own spelling, and an exact comparison silently dropped the mapping (the
// temperature then stayed at the SOLUTION value).
static bool phreeqc_same_name(const std::string& a, const std::string& b){
    return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(),
        [](char x, char y){ return std::toupper(static_cast<unsigned char>(x))
                                == std::toupper(static_cast<unsigned char>(y)); });
}

static std::string phreeqc_resolve_path(const std::string& path) {
    if (!path.empty() && path[0] == '/') return path;
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) return std::string(cwd) + "/" + path;
    return path;
}

// Value of a CALCULATE_VALUES function of the .pqi file (RENEWAL time scales
// given by name, so that they are calibrated like any other .pqi constant).
// Evaluated in a separate IPhreeqc instance: nothing changes in PhreeqcRM.
static bool phreeqc_calc_value(const std::string& database, const std::string& pqi,
                               const std::string& name, double& value) {
    IPhreeqc ip;
    if (ip.LoadDatabase(database.c_str()) != 0) return false;
    if (ip.RunFile(pqi.c_str()) != 0) return false;
    const std::string in =
        "SELECTED_OUTPUT 99999\n    -reset false\n"
        "USER_PUNCH 99999\n    -headings v\n    10 PUNCH CALC_VALUE(\"" + name + "\")\n"
        "SOLUTION 99999\nEND\n";
    if (ip.RunString(in.c_str()) != 0) return false;
    ip.SetCurrentSelectedOutputUserNumber(99999);
    if (ip.GetSelectedOutputRowCount() < 2) return false;
    VAR v; VarInit(&v);
    ip.GetSelectedOutputValue(1, 0, &v);
    const bool ok = (v.type == TT_DOUBLE);
    if (ok) value = v.dVal;
    VarClear(&v);
    return ok;
}

/* #################################################
// One-time preparation of the coupling (first call of phreeqc_run):
// (1) the carrier species of every mobile component (charge, H and O that
//     move with a transported element total), and
// (2) the RENEWAL settings of each compartment: reference solution and time
//     scale per cell, and
// (3) the read-out components (DIAGNOSTICS).
################################################# */
static void phreeqc_prepare_coupling(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output){

    auto& P = *OpenWQ_wqconfig.CH_model->PHREEQC;
    PhreeqcRM& rm = *P.phreeqcrm;
    const int nxyz = rm.GetGridCellCount();
    const std::vector<std::string>& comps = rm.GetComponents();
    std::string msg_string;

    for (size_t j = 0; j < comps.size(); j++) {
        if (comps[j] == "H") P.comp_H = (int)j;
        else if (comps[j] == "O") P.comp_O = (int)j;
        else if (comps[j] == "Charge") P.comp_Charge = (int)j;
        else if (comps[j] == "H2O") P.comp_H2O = (int)j;
    }
    if (P.comp_H < 0 || P.comp_O < 0 || P.comp_Charge < 0) {
        msg_string = "<OpenWQ> ERROR: PHREEQC components H, O and Charge not found.";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        exit(EXIT_FAILURE);
    }

    // (1) carriers
    const std::vector<std::string>& sp_names = rm.GetSpeciesNames();
    const std::vector<cxxNameDouble>& sp_stoich = rm.GetSpeciesStoichiometry();
    P.carriers.assign(P.mobile_species.size(), {});
    for (size_t m = 0; m < P.mobile_species.size(); m++) {
        const std::string& elt = comps[P.mobile_species[m]];
        std::decay_t<decltype(P.carriers)>::value_type simple, all;
        for (size_t s = 0; s < sp_stoich.size(); s++) {
            double nu = 0.0, h = 0.0, o = 0.0, z = 0.0;
            bool other = false;
            for (const auto& e : sp_stoich[s]) {
                if (e.first == elt) nu = e.second;
                else if (e.first == "H") h = e.second;
                else if (e.first == "O") o = e.second;
                else if (e.first == "Charge") z = e.second;
                else other = true;
            }
            if (nu <= 0.0) continue;
            all.push_back({(int)s, nu, h, o, z});
            if (!other) simple.push_back({(int)s, nu, h, o, z});
        }
        P.carriers[m] = simple.empty() ? all : simple;
        // Optional carriers named in the module JSON, "CARRIERS": {"C": ["CO2"]}.
        // The default weights the carriers by the species of the destination
        // cell, which is exact for a component whose charge does not depend on
        // the speciation (nitrate). Inorganic carbon arrives from the soil as
        // CO2 but is weighted as HCO3- in a stream at pH 8, which adds a charge
        // that the water did not bring and raises the pH. Carried as CO2, the
        // carbon moves without charge, the alkalinity moves with the charged
        // components (Ca+2), and each cell speciates its carbon from its own
        // charge balance (transport of DIC and alkalinity).
        if (OpenWQ_json.BGC_module.contains("CARRIERS")) {
            for (const auto& it : OpenWQ_json.BGC_module["CARRIERS"].items()) {
                if (!phreeqc_same_name(it.key(), elt)) continue;
                std::decay_t<decltype(P.carriers)>::value_type named;
                for (const auto& v : it.value())
                    for (const auto& cs : all)
                        if (phreeqc_same_name(sp_names[cs.s], v.get<std::string>())) named.push_back(cs);
                if (named.empty()) {
                    msg_string = "<OpenWQ> ERROR: PHREEQC CARRIERS of " + elt + ": none of the named species contains " + elt;
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    exit(EXIT_FAILURE);
                }
                P.carriers[m] = named;
            }
        }
        msg_string = "<OpenWQ> PHREEQC: mobile component " + elt + " carries";
        if (simple.empty()) msg_string += " (no simple species, all its species)";
        for (const auto& cs : P.carriers[m])
            msg_string += " " + sp_names[cs.s] + " (H " + std::to_string(cs.h)
                + ", O " + std::to_string(cs.o) + ", charge " + std::to_string(cs.z)
                + " per " + std::to_string(cs.nu) + ")";
        if (P.carriers[m].empty()) msg_string += " nothing (no species found: written without charge, H, O)";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

    P.neutral_carriers = OpenWQ_json.BGC_module.contains("NEUTRAL_CARRIERS")
        && OpenWQ_json.BGC_module["NEUTRAL_CARRIERS"].get<bool>();
    if (P.neutral_carriers) {
        msg_string = "<OpenWQ> PHREEQC: NEUTRAL_CARRIERS: every mobile component moves as the neutral formula of its "
            "species (charge balanced with OH- or H+); exact transport of electroneutral water when all charged "
            "components are mobile";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

    // RESET_WATER: the water of a cell is the host's and is not transported, so the water that the reactions and the
    // carriers produce in a cell stays there when nothing renews the cell (a reach that loses CO2 while the transport
    // refills its carbon as H2CO3 gains one H2O per mol). With "RESET_WATER": true the H2O component of every cell is
    // set back to its value at the first call before each RunCells; the excess H, O and charge are kept.
    {
        const auto& bm = OpenWQ_json.BGC_module;
        P.reset_water = bm.contains("RESET_WATER") && bm["RESET_WATER"].is_boolean() && bm["RESET_WATER"].get<bool>();
    }
    if (P.reset_water) {
        if (P.comp_H2O < 0) {
            msg_string = "<OpenWQ> ERROR: PHREEQC RESET_WATER needs the H2O component (COMPONENT_H2O true).";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
            exit(EXIT_FAILURE);
        }
        msg_string = "<OpenWQ> PHREEQC: RESET_WATER: the water (H2O component) of every cell is set to its initial value "
            "before each RunCells (the water of the host volume)";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

    // (2) renewal
    P.renewal_tau_s.assign(nxyz, 0.0);
    P.renewal_ref.assign(nxyz, -1);
    P.renewal_ref2.assign(nxyz, -1);
    P.renewal_mobile.assign(nxyz, 0);
    P.renewal_tracer.assign(nxyz, -1);
    P.renewal_g0.assign(nxyz, 0.0);
    P.renewal_g1.assign(nxyz, 1.0);
    P.renewal_age_tracer.assign(nxyz, -1);
    P.renewal_age_unit_s.assign(nxyz, 0.0);
    P.renewal_age_min_s.assign(nxyz, 0.0);
    P.renewal_ref_c.clear();
    std::map<int, int> ref_of_solution;
    const std::string database = phreeqc_resolve_path(OpenWQ_json.BGC_module["DATABASE"].get<std::string>());
    const std::string pqi = phreeqc_resolve_path(OpenWQ_json.BGC_module["FILEPATH"].get<std::string>());
    int cell0 = 0;
    for (unsigned int icmp = 0; icmp < OpenWQ_hostModelconfig.get_num_HydroComp(); icmp++) {
        const int ncell = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp)
                        * OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp)
                        * OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);
        const std::string cname = OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);
        const json& bgc = OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"];
        if (bgc.contains(cname) && bgc[cname].contains("RENEWAL")) {
            const json& ren = bgc[cname]["RENEWAL"];
            int sol = -1;
            double tau = 0.0;
            std::string tau_label;
            // TIME_SCALE {"TRACER": <component>, "UNIT_S": s, "MIN_S": s}: the time scale of
            // each cell is the concentration of a transported age tracer (mg/L) x UNIT_S,
            // at least MIN_S. With a tracer produced at one unit per UNIT_S in the reaches
            // and absent from their inflows (the time the water has spent in the river
            // network), the renewal is the replacement of the reach water by its inflow at
            // the rate of a well-mixed reach, 1 / residence time: young water at high flow
            // stays close to the composition of its inflow, old water at low flow has had
            // the time to exchange CO2 with the atmosphere and to photosynthesize.
            int age_tracer = -1;
            double age_unit = 0.0, age_min = 0.0;
            try {
                sol = ren.at("SOLUTION").get<int>();
                const json& ts = ren.at("TIME_SCALE");
                if (ts.is_object()) {
                    const std::string aname = ts.at("TRACER").get<std::string>();
                    age_unit = ts.at("UNIT_S").get<double>();
                    age_min = ts.contains("MIN_S") ? ts["MIN_S"].get<double>() : 60.0;
                    for (size_t j = 0; j < comps.size(); j++)
                        if (phreeqc_same_name(comps[j], aname)) age_tracer = (int)j;
                    if (age_tracer < 0 || !(age_unit > 0.0) || !(age_min > 0.0)) {
                        msg_string = "<OpenWQ> ERROR: RENEWAL of " + cname + ": TIME_SCALE TRACER " + aname
                            + " is not a component, or UNIT_S / MIN_S are not > 0";
                        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                        exit(EXIT_FAILURE);
                    }
                    tau = age_min;
                    tau_label = "age tracer " + comps[age_tracer] + " x " + std::to_string(age_unit) + " s per mg/L (at least "
                        + std::to_string(age_min) + " s)";
                } else if (ts.is_number()) {
                    tau = ts.get<double>();
                    tau_label = std::to_string(tau) + " s";
                } else {
                    const std::string name = ts.get<std::string>();
                    if (!phreeqc_calc_value(database, pqi, name, tau)) {
                        msg_string = "<OpenWQ> ERROR: RENEWAL of " + cname + ": CALCULATE_VALUES '"
                            + name + "' could not be evaluated from " + pqi;
                        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                        exit(EXIT_FAILURE);
                    }
                    tau_label = name + " = " + std::to_string(tau) + " s";
                }
            } catch (...) {
                msg_string = "<OpenWQ> ERROR: RENEWAL of " + cname
                    + " needs {\"SOLUTION\": <id>, \"TIME_SCALE\": <seconds or CALCULATE_VALUES name>}";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                exit(EXIT_FAILURE);
            }
            if (!(tau > 0.0) || !std::isfinite(tau)) {
                msg_string = "<OpenWQ> ERROR: RENEWAL of " + cname + ": time scale must be > 0 (" + tau_label + ")";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                exit(EXIT_FAILURE);
            }
            if (!ref_of_solution.count(sol)) {
                std::vector<double> cref;
                std::vector<int> bc1(1, sol);
                if (rm.InitialPhreeqc2Concentrations(cref, bc1) != IRM_OK
                        || cref.size() != comps.size()) {
                    msg_string = "<OpenWQ> ERROR: RENEWAL of " + cname + ": SOLUTION "
                        + std::to_string(sol) + " not found in " + pqi;
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    exit(EXIT_FAILURE);
                }
                ref_of_solution[sol] = (int)P.renewal_ref_c.size();
                P.renewal_ref_c.push_back(cref);
            }
            // Optional end-member mixing, "MIX": {"SOLUTION": <id2>, "TRACER": <component>}:
            // the reference becomes (1 - g) SOLUTION + g SOLUTION id2, with g the
            // concentration of the tracer (a transported conservative component,
            // e.g. a ground-water tracer held at its value of SOLUTION id2 in the
            // aquifer) relative to its value in SOLUTION id2, bounded to [0, 1].
            // A mixture of two solutions in component space is a solution (PHREEQC
            // MIX), so that the reach chemistry follows the mixing of the waters
            // that the host model routes to it.
            int ref2 = -1, tracer = -1;
            double g0 = 0.0, g1 = 1.0;
            std::string mix_label;
            if (ren.contains("MIX")) {
                int sol2 = -1;
                std::string tname;
                try {
                    sol2 = ren["MIX"].at("SOLUTION").get<int>();
                    tname = ren["MIX"].at("TRACER").get<std::string>();
                } catch (...) {
                    msg_string = "<OpenWQ> ERROR: RENEWAL of " + cname
                        + ": MIX needs {\"SOLUTION\": <id>, \"TRACER\": <component>}";
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    exit(EXIT_FAILURE);
                }
                for (size_t j = 0; j < comps.size(); j++)
                    if (phreeqc_same_name(comps[j], tname)) tracer = (int)j;
                if (!ref_of_solution.count(sol2)) {
                    std::vector<double> cref2;
                    std::vector<int> bc2(1, sol2);
                    if (rm.InitialPhreeqc2Concentrations(cref2, bc2) != IRM_OK
                            || cref2.size() != comps.size()) {
                        msg_string = "<OpenWQ> ERROR: RENEWAL of " + cname + ": MIX SOLUTION "
                            + std::to_string(sol2) + " not found in " + pqi;
                        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                        exit(EXIT_FAILURE);
                    }
                    ref_of_solution[sol2] = (int)P.renewal_ref_c.size();
                    P.renewal_ref_c.push_back(cref2);
                }
                ref2 = ref_of_solution[sol2];
                if (tracer < 0 || !(P.renewal_ref_c[ref2][tracer] > 0.0)) {
                    msg_string = "<OpenWQ> ERROR: RENEWAL of " + cname + ": MIX TRACER " + tname
                        + " is not a component or is zero in SOLUTION " + std::to_string(sol2);
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    exit(EXIT_FAILURE);
                }
                // Optional "RANGE": [g0, g1]: tracer fractions taken as the first and
                // the second end-member (the tracer of a reach never reaches 0 or 1
                // when it is a contact-time tracer of mixed waters)
                if (ren["MIX"].contains("RANGE")) {
                    try {
                        g0 = ren["MIX"]["RANGE"].at(0).get<double>();
                        g1 = ren["MIX"]["RANGE"].at(1).get<double>();
                    } catch (...) { g0 = 0.0; g1 = 0.0; }
                    if (!(g1 > g0)) {
                        msg_string = "<OpenWQ> ERROR: RENEWAL of " + cname + ": MIX RANGE needs [g0, g1] with g1 > g0";
                        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                        exit(EXIT_FAILURE);
                    }
                }
                mix_label = ", mixed with SOLUTION " + std::to_string(sol2) + " by the fraction of " + comps[tracer]
                    + " (fractions " + std::to_string(g0) + " to " + std::to_string(g1) + " of its value there)";
            }
            for (int i = 0; i < ncell; i++) {
                P.renewal_tau_s[cell0 + i] = tau;
                P.renewal_ref[cell0 + i] = ref_of_solution[sol];
                P.renewal_ref2[cell0 + i] = ref2;
                P.renewal_tracer[cell0 + i] = tracer;
                P.renewal_g0[cell0 + i] = g0;
                P.renewal_g1[cell0 + i] = g1;
                P.renewal_age_tracer[cell0 + i] = age_tracer;
                P.renewal_age_unit_s[cell0 + i] = age_unit;
                P.renewal_age_min_s[cell0 + i] = age_min;
                P.renewal_mobile[cell0 + i] = (ren.contains("MOBILE") && ren["MOBILE"].get<bool>()) ? 1 : 0;
            }
            msg_string = "<OpenWQ> PHREEQC: RENEWAL of " + cname + " towards SOLUTION "
                + std::to_string(sol) + mix_label + ", time scale " + tau_label
                + ((ren.contains("MOBILE") && ren["MOBILE"].get<bool>()) ? " (MOBILE: transported components too)" : "");
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }
        cell0 += ncell;
    }

    // (3) read-out components (DIAGNOSTICS)
    P.diagnostics.clear();
    P.is_diagnostic.assign(comps.size(), 0);
    P.species_Hplus = -1;
    for (size_t s = 0; s < sp_names.size(); s++)
        if (sp_names[s] == "H+") P.species_Hplus = (int)s;
    if (OpenWQ_json.BGC_module.contains("DIAGNOSTICS")) {
        const json& dg = OpenWQ_json.BGC_module["DIAGNOSTICS"];
        for (auto it = dg.begin(); it != dg.end(); ++it) {
            std::string what;
            try { what = it.value().get<std::string>(); } catch (...) { what.clear(); }
            int jc = -1, mob = -1, sp = -1;
            for (size_t j = 0; j < comps.size(); j++)
                if (phreeqc_same_name(comps[j], it.key())) jc = (int)j;
            for (size_t m = 0; m < P.mobile_species.size(); m++)
                if ((int)P.mobile_species[m] == jc) mob = (int)m;
            const bool is_ph = phreeqc_same_name(what, "PH");
            if (!is_ph) {
                // JSON values arrive upper-cased: species matched without regard
                // to case, an exact match preferred
                int n_match = 0;
                for (size_t s = 0; s < sp_names.size(); s++) {
                    if (sp_names[s] == what) { sp = (int)s; n_match = 1; break; }
                    if (phreeqc_same_name(sp_names[s], what)) { sp = (int)s; n_match++; }
                }
                if (n_match > 1) sp = -2;
            }
            std::string err;
            if (jc < 0) err = "component '" + it.key() + "' not found (a read-out component must be"
                " nonzero in a SOLUTION of the .pqi file to be a PhreeqcRM component)";
            else if (mob < 0) err = "component '" + it.key() + "' is not in BGC_GENERAL_MOBILE_SPECIES (no output)";
            else if (is_ph && P.species_Hplus < 0) err = "species H+ not found";
            else if (!is_ph && sp == -1) err = "aqueous species '" + what + "' not found";
            else if (!is_ph && sp == -2) err = "aqueous species '" + what + "' is ambiguous (names differ only by case)";
            if (!err.empty()) {
                msg_string = "<OpenWQ> ERROR: PHREEQC DIAGNOSTICS: " + err
                    + ". Use {\"<component>\": \"pH\" or \"<aqueous species>\"}.";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                exit(EXIT_FAILURE);
            }
            P.diagnostics.push_back({jc, is_ph ? -1 : sp});
            P.is_diagnostic[jc] = 1;
            msg_string = "<OpenWQ> PHREEQC: read-out component " + comps[jc] + " = "
                + (is_ph ? std::string("pH (-log10 a H+)") : "molality of " + sp_names[sp])
                + " after each RunCells (not part of the chemistry)";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }
    }
    // (4) host inputs (HOST_INPUTS): components set from host variables before every RunCells
    P.host_inputs.clear();
    if (OpenWQ_json.BGC_module.contains("HOST_INPUTS")) {
        const json& hi = OpenWQ_json.BGC_module["HOST_INPUTS"];
        for (auto it = hi.begin(); it != hi.end(); ++it) {
            int jc = -1;
            for (size_t j = 0; j < comps.size(); j++)
                if (phreeqc_same_name(comps[j], it.key())) jc = (int)j;
            // flat form, as DIAGNOSTICS: {"<component>": "RELATIVE_VOLUME"}; the time scale of the running mean is
            // the top-level "HOST_INPUTS_TAU_DAYS" (365 days by default)
            std::string var;
            double tau_d = 365.0;
            try { var = it.value().get<std::string>(); } catch (...) { var.clear(); }
            try {
                if (OpenWQ_json.BGC_module.contains("HOST_INPUTS_TAU_DAYS"))
                    tau_d = OpenWQ_json.BGC_module["HOST_INPUTS_TAU_DAYS"].get<double>();
            } catch (...) { tau_d = 365.0; }
            if (jc < 0 || !phreeqc_same_name(var, "RELATIVE_VOLUME") || !(tau_d > 0.0)) {
                msg_string = "<OpenWQ> ERROR: PHREEQC HOST_INPUTS " + it.key() + ": the component must be nonzero in a SOLUTION"
                    " of the .pqi file, and the value RELATIVE_VOLUME ({\"<component>\": \"RELATIVE_VOLUME\"})";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                exit(EXIT_FAILURE);
            }
            P.host_inputs.push_back({jc, tau_d * 86400.0});
            msg_string = "<OpenWQ> PHREEQC: host input component " + comps[jc] + " = water volume / its running geometric mean ("
                + std::to_string(tau_d) + " days), set before each RunCells";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }
    }
    P.host_vol_mean.assign(nxyz, 0.0);
    P.coupling_ready = true;
}

/* #################################################
// Compute each chemical transformation using PHREEQC
//
// UNIT CONVERSION NOTES:
// ----------------------
// OpenWQ internal units:
//   - chemass: stored in grams (g)
//   - water volume of the host: m3, so chemass/volume is g/m3 = mg/L
//   - volumes passed to the conversions below: litres (L) = host m3 * 1000
//
// PHREEQC PhreeqcRM units (with SetUnitsSolution(2)):
//   - concentration: mol/kgw (moles per kilogram water)
//   - For dilute solutions: mol/kgw ≈ mol/L (since water density ~1 kg/L)
//
// Conversion formulas:
//   mg/L to mol/kgw:  conc_mol = conc_mg_L / (GFW * 1000)
//   mol/kgw to mg/L:  conc_mg_L = conc_mol * GFW * 1000
//
// where GFW = Gram Formula Weight in g/mol
//
################################################# */
void OpenWQ_CH_model::phreeqc_run(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output){

    std::vector<double> c, volumes, temperatures, pressures;
    int nxyz;
    int nx, ny, nz;
    int indx;
    std::string dependancy;
    std::string compName_icmp;
    unsigned int ndep;
    int temperature_indx;
    int pressure_indx;
    std::string msg_string;

    nxyz = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetGridCellCount();

    // Get the actual component count from PhreeqcRM
    int nc = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetComponentCount();
    msg_string = "<OpenWQ> PHREEQC phreeqc_run: nxyz=" + std::to_string(nxyz)
        + ", nc=" + std::to_string(nc)
        + ", num_chem=" + std::to_string(OpenWQ_wqconfig.CH_model->PHREEQC->num_chem);
    if (OpenWQ_wqconfig.debug_mode) OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    // Get GFW for unit conversion
    const std::vector<double>& gfw = OpenWQ_wqconfig.CH_model->PHREEQC->gfw;
    bool has_gfw = (gfw.size() >= (size_t)nc);

    if (!has_gfw) {
        msg_string = "<OpenWQ> WARNING: GFW not available for all components. "
            "Unit conversion disabled (assuming units already in mol/kgw).";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

    // Use nc (actual PhreeqcRM component count) for sizing
    c.resize(nxyz * nc);

    // Pre-load the CURRENT PHREEQC concentrations for ALL components before we
    // overwrite the mobile species below. Without this, every NON-mobile
    // component (H, Charge, and any ion openWQ does not transport) would be
    // pushed to PHREEQC as 0 by SetConcentrations, destroying the water and
    // charge balance and guaranteeing the geochemical solve fails to converge.
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetConcentrations(c);

    auto& P = *OpenWQ_wqconfig.CH_model->PHREEQC;
    if (!P.coupling_ready)
        phreeqc_prepare_coupling(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_output);
    const auto& cell_is_active = P.cell_active;
    auto active_cell = [&cell_is_active](int i){ return cell_is_active.empty() || cell_is_active[i] != 0; };

    // Charge, H and O per mole of each mobile component in each cell, from the
    // species after the previous RunCells (the state that c holds now). A cell
    // without the component takes the domain-wide speciation, and a domain
    // without it the first carrier species.
    const unsigned int n_mob = P.mobile_species.size();
    std::vector<double> carry_h(n_mob * nxyz, 0.0), carry_o(n_mob * nxyz, 0.0), carry_z(n_mob * nxyz, 0.0);
    {
        std::vector<double> sc;
        P.phreeqcrm->GetSpeciesConcentrations(sc);
        const size_t nsp = P.phreeqcrm->GetSpeciesCount();
        const bool have_sc = (sc.size() == nsp * (size_t)nxyz);
        for (unsigned int m = 0; m < n_mob && m < P.carriers.size(); m++) {
            const auto& cs = P.carriers[m];
            if (cs.empty()) continue;
            double dom_k = 0.0, dom_h = 0.0, dom_o = 0.0, dom_z = 0.0;
            std::vector<double> cell_k(nxyz, 0.0);
            for (int i = 0; i < nxyz; i++) {
                if (!active_cell(i) || !have_sc) continue;
                double sk = 0.0, sh = 0.0, so = 0.0, sz = 0.0;
                for (const auto& x : cs) {
                    const double ms = std::max(0.0, sc[(size_t)x.s * nxyz + i]);
                    sk += x.nu * ms; sh += x.h * ms; so += x.o * ms; sz += x.z * ms;
                }
                cell_k[i] = sk;
                if (sk > 0.0) {
                    carry_h[m * nxyz + i] = sh / sk; carry_o[m * nxyz + i] = so / sk; carry_z[m * nxyz + i] = sz / sk;
                }
                dom_k += sk; dom_h += sh; dom_o += so; dom_z += sz;
            }
            double def_h, def_o, def_z;
            if (dom_k > 0.0) { def_h = dom_h / dom_k; def_o = dom_o / dom_k; def_z = dom_z / dom_k; }
            else { def_h = cs[0].h / cs[0].nu; def_o = cs[0].o / cs[0].nu; def_z = cs[0].z / cs[0].nu; }
            for (int i = 0; i < nxyz; i++) {
                if (cell_k[i] > 0.0) continue;
                carry_h[m * nxyz + i] = def_h; carry_o[m * nxyz + i] = def_o; carry_z[m * nxyz + i] = def_z;
            }
        }
    }
    // NEUTRAL_CARRIERS: the charge of the carrier balanced with OH- (cations: H and O + z) or H+ (anions: H + |z|).
    // Every species of an element then carries the same formula up to water, whatever the speciation of the cell.
    if (P.neutral_carriers) {
        for (size_t q = 0; q < carry_z.size(); q++) {
            const double z = carry_z[q];
            carry_h[q] += std::fabs(z);
            carry_o[q] += std::max(z, 0.0);
            carry_z[q] = 0.0;
        }
    }

    // Renewal of the non-transported chemistry (all components; the mobile
    // ones are overwritten below, carrying their charge, H and O, so that the
    // result equals renewing the non-transported part only).
    if (OpenWQ_hostModelconfig.get_time_step() > 0 && !P.renewal_ref_c.empty()) {
        const double dt = OpenWQ_hostModelconfig.get_time_step();
        const std::vector<double>& gfw_r = OpenWQ_wqconfig.CH_model->PHREEQC->gfw;
        const bool has_gfw_r = (gfw_r.size() >= (size_t)nc);
        for (int i = 0; i < nxyz; i++) {
            if (!active_cell(i) || P.renewal_ref[i] < 0) continue;
            double tau_i = P.renewal_tau_s[i];
            if (P.renewal_age_tracer[i] >= 0) {
                // time scale from the age tracer of the cell (mol/L -> mg/L with its GFW)
                const int ka = P.renewal_age_tracer[i];
                const double gfw_a = (has_gfw_r && gfw_r[ka] > 0.0) ? gfw_r[ka] : 1.0;
                tau_i = std::max(P.renewal_age_min_s[i], c[(size_t)ka * nxyz + i] * gfw_a * 1000.0 * P.renewal_age_unit_s[i]);
            }
            const double f = 1.0 - std::exp(-dt / tau_i);
            const std::vector<double>& cref = P.renewal_ref_c[P.renewal_ref[i]];
            if (P.renewal_ref2[i] >= 0) {
                // MIX: fraction of the second end-member from the tracer of the
                // cell (its transported value of the previous step)
                const std::vector<double>& cref2 = P.renewal_ref_c[P.renewal_ref2[i]];
                const int k = P.renewal_tracer[i];
                const double g = std::min(1.0, std::max(0.0,
                    (c[(size_t)k * nxyz + i] / cref2[k] - P.renewal_g0[i]) / (P.renewal_g1[i] - P.renewal_g0[i])));
                for (int j = 0; j < nc; j++)
                    c[(size_t)j * nxyz + i] += f * ((1.0 - g) * cref[j] + g * cref2[j] - c[(size_t)j * nxyz + i]);
                continue;
            }
            for (int j = 0; j < nc; j++)
                c[(size_t)j * nxyz + i] += f * (cref[j] - c[(size_t)j * nxyz + i]);
        }
    }

    volumes.resize(nxyz);
    temperatures.resize(nxyz);
    pressures.resize(nxyz);
    indx = 0;
    pressures = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetPressure();
    temperatures = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetTemperature();
    // Lowest temperature passed to PHREEQC in each cell (C): 0 C unless the
    // compartment's TEMPERATURE entry sets "MIN_C" (see the clamp below).
    std::vector<double> temperatures_min(nxyz, 0.0);
    if (P.temp_tau_s.size() != (size_t)nxyz) { P.temp_tau_s.assign(nxyz, 0.0); P.temp_mean.assign(nxyz, std::nan("")); }

    // First pass: collect volumes, temperatures, pressures
    for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){


        nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
        ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
        nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements
        compName_icmp = OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);
        // Water temperature of the compartment: a host variable ("Tsoil_K"),
        // or a linear relation to one, {"DEPENDENCY": "Tair_K", "LINEAR": [a, b]}
        // for T = a + b * variable (C), e.g. stream temperature from air
        // temperature. Below 0 C it is set to 0 C further down.
        double temp_a = 0.0, temp_b = 1.0, temp_min = 0.0, temp_tau = 0.0;
        try {
            const json& temp_map = OpenWQ_json.BGC_module["TEMPERATURE"][compName_icmp];
            if (temp_map.is_object()) {
                dependancy = temp_map.at("DEPENDENCY").get<std::string>();
                if (temp_map.contains("LINEAR")) {
                    temp_a = temp_map["LINEAR"].at(0).get<double>();
                    temp_b = temp_map["LINEAR"].at(1).get<double>();
                }
                if (temp_map.contains("MIN_C"))
                    temp_min = temp_map["MIN_C"].get<double>();
                if (temp_map.contains("SMOOTH_HOURS"))
                    temp_tau = 3600.0 * temp_map["SMOOTH_HOURS"].get<double>();
            } else {
                dependancy = temp_map.get<std::string>();
            }
            temperature_indx = -1;
            for (ndep = 0; ndep < OpenWQ_hostModelconfig.get_num_HydroDepend(); ndep ++) {
                if (phreeqc_same_name(OpenWQ_hostModelconfig.get_HydroDepend_name_at(ndep), dependancy)) temperature_indx = ndep;
            }
            if (temperature_indx == -1) {
                // Create Message
                msg_string = "<OpenWQ> Unknown TEMPERATURE dependancy with name '"
                    + dependancy
                    + "' defined for compartment: "
                    + compName_icmp;

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(
                    OpenWQ_wqconfig,    // for Log file name
                    msg_string,         // message
                    true,               // print in console
                    true);              // print in log file

            }
        } catch(...) {
            temperature_indx = -1;
        }

        // Pressure of the compartment (atm): a host variable, or a linear relation
        // to one, {"DEPENDENCY": "SWrad_Wm2", "LINEAR": [a, b]} for P = a + b * variable.
        // With a = 1 and a small b it carries a host variable that PHREEQC has no
        // other access to (e.g. light for in-stream photosynthesis, read in RATES as
        // (PRESSURE - a) / b) while changing the equilibrium constants negligibly.
        double press_a = 0.0, press_b = 1.0;
        try {
            const json& press_map = OpenWQ_json.BGC_module["PRESSURE"][compName_icmp];
            if (press_map.is_object()) {
                dependancy = press_map.at("DEPENDENCY").get<std::string>();
                if (press_map.contains("LINEAR")) {
                    press_a = press_map["LINEAR"].at(0).get<double>();
                    press_b = press_map["LINEAR"].at(1).get<double>();
                }
            } else {
                dependancy = press_map.get<std::string>();
            }
            pressure_indx = -1;
            for (ndep = 0; ndep < OpenWQ_hostModelconfig.get_num_HydroDepend(); ndep ++) {
                if (phreeqc_same_name(OpenWQ_hostModelconfig.get_HydroDepend_name_at(ndep), dependancy)) pressure_indx = ndep;
            }
            if (pressure_indx == -1) {
                // Create Message
                msg_string = "<OpenWQ> Unknown PRESSURE dependancy with name '"
                    + dependancy
                    + "' defined for compartment: "
                    + compName_icmp;

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(
                    OpenWQ_wqconfig,    // for Log file name
                    msg_string,         // message
                    true,               // print in console
                    true);              // print in log file

            }
        } catch(...) {
            pressure_indx = -1;
        }
        for (int ix=0; ix<nx; ix++) {
            for (int iy = 0; iy<ny;iy++) {
                for (int iz=0; iz<nz; iz++) {
                    // Host water volume is in m3; the conversions below use
                    // litres. Treating m3 as litres made every concentration
                    // handed to PHREEQC 1000 times too high (the mass round trip
                    // hid it, but equilibria and rate laws saw the wrong values).
                    volumes[indx] = 1000.0 * OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(icmp, ix, iy, iz);
                    if (volumes[indx] == 0) volumes[indx] = 1;

                    if (temperature_indx!=-1) {
                        double temp_val = OpenWQ_hostModelconfig.get_dependVar_at(temperature_indx, ix, iy, iz);
                        // Auto-detect if temperature is in Kelvin (>100 suggests Kelvin)
                        // PHREEQC expects Celsius, so convert if needed
                        // Typical water temps: 0-40°C or 273-313K
                        if (temp_val > 100.0) {
                            // Likely Kelvin - convert to Celsius
                            temp_val = temp_val - 273.15;
                        }
                        temperatures[indx] = temp_a + temp_b * temp_val;
                    }
                    temperatures_min[indx] = temp_min;
                    P.temp_tau_s[indx] = (temperature_indx != -1) ? temp_tau : 0.0;
                    if (pressure_indx != -1) {
                        pressures[indx] = press_a + press_b * OpenWQ_hostModelconfig.get_dependVar_at(pressure_indx, ix, iy, iz);
                    }
                    indx++;

                }
            }
        }
    }


    // ########################################################
    // UNIT CONVERSION: OpenWQ (mg/L) -> PHREEQC (mol/kgw)
    // ########################################################
    // OpenWQ chemass is in grams (g)
    // concentration = chemass / volume gives g/L
    // To get mg/L: multiply by 1000
    // To get mol/kgw: divide mg/L by (GFW * 1000)
    // Combined: conc_mol = (chemass_g / volume_L) / GFW
    //         = chemass_g / (volume_L * GFW)
    // ########################################################

    // Use number of mobile species from OpenWQ
    unsigned int num_mobile = OpenWQ_wqconfig.CH_model->PHREEQC->mobile_species.size();
    unsigned int num_chem_to_use = std::min(num_mobile, OpenWQ_wqconfig.CH_model->PHREEQC->num_chem);

    for (unsigned int chemi=0;chemi<num_chem_to_use;chemi++){
        indx = 0;

        // Get the PHREEQC component index for this OpenWQ chemical
        unsigned int phreeqc_idx = OpenWQ_wqconfig.CH_model->PHREEQC->mobile_species[chemi];

        // Get GFW for this PHREEQC component (default to 1.0 if not available)
        double component_gfw = (has_gfw && phreeqc_idx < gfw.size()) ? gfw[phreeqc_idx] : 1.0;

        // Protect against invalid GFW
        if (component_gfw <= 0 || std::isnan(component_gfw) || std::isinf(component_gfw)) {
            component_gfw = 1.0;
        }

        for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){

            nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
            ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
            nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements
            for (int ix=0; ix<nx; ix++) {
                for (int iy = 0; iy<ny;iy++) {
                    for (int iz=0; iz<nz; iz++) {
                        // OpenWQ: chemass is in grams, volume in liters
                        // conc_g_L = chemass_g / volume_L
                        // Index chemass by the PHREEQC component index (phreeqc_idx),
                        // NOT the mobile ordinal (chemi) — this matches how transport
                        // (native_adv/advdisp), IC and output index chemass. Using chemi
                        // here read the wrong slot, feeding ~0 to PHREEQC.
                        double conc_g_L = (*OpenWQ_vars.chemass)(icmp)(phreeqc_idx)(ix,iy,iz) / volumes[indx];

                        // Convert g/L to mol/kgw (mol/L for dilute solutions)
                        // mol/L = g/L / GFW(g/mol)
                        double conc_mol_kgw = conc_g_L / component_gfw;

                        // Validate - replace NaN/Inf with 0
                        if (std::isnan(conc_mol_kgw) || std::isinf(conc_mol_kgw) || conc_mol_kgw < 0) {
                            conc_mol_kgw = 0.0;
                        }

                        // A read-out component is not part of the chemistry
                        // (set after RunCells, see DIAGNOSTICS)
                        const bool readout = P.is_diagnostic[phreeqc_idx] != 0;
                        if (readout) conc_mol_kgw = 0.0;

                        // Use PHREEQC component index for concentration array.
                        // The change moves the charge, H and O of the
                        // component's species with it (see carriers).
                        if (active_cell(indx) && !readout) {
                            const double dk = conc_mol_kgw - c[(size_t)phreeqc_idx*nxyz+indx];
                            c[(size_t)P.comp_H*nxyz+indx]      += carry_h[chemi*nxyz+indx] * dk;
                            c[(size_t)P.comp_O*nxyz+indx]      += carry_o[chemi*nxyz+indx] * dk;
                            c[(size_t)P.comp_Charge*nxyz+indx] += carry_z[chemi*nxyz+indx] * dk;
                        }
                        c[phreeqc_idx*nxyz+indx] = conc_mol_kgw;
                        indx++;
                    }
                }
            }
        }
    }

    // HOST_INPUTS: water volume of the cell / its running mean (mg/L of the component = the ratio, GFW 1)
    if (!P.host_inputs.empty()) {
        const double dt = OpenWQ_hostModelconfig.get_time_step();
        for (const auto& hin : P.host_inputs) {
            const double w = (dt > 0) ? 1.0 - std::exp(-dt / hin.tau_s) : 0.0;
            for (int i = 0; i < nxyz; i++) {
                if (!active_cell(i)) continue;
                const double v = volumes[i];
                if (!(v > 0.0) || !std::isfinite(v)) continue;
                // running GEOMETRIC mean (the typical volume): the volume of a reach is skewed by floods, and an
                // arithmetic mean sits far above the volume of most days
                double& m = P.host_vol_mean[i];
                if (!(m > 0.0)) m = v; else m = std::exp(std::log(m) + w * (std::log(v) - std::log(m)));
                {
                    static const bool dbg = std::getenv("OPENWQ_DEBUG_VOLR") != nullptr;
                    const long step = OpenWQ_hostModelconfig.get_interaction_step();
                    if (dbg && (step <= 24 || step % 2000 == 0))
                        std::fprintf(stderr, "VOLR phreeqc %ld %.10e %.10e %.10e\n", step, v / 1000.0, m / 1000.0, v / m);
                }
                c[(size_t)hin.comp * nxyz + i] = (v / m) / 1000.0;
            }
        }
    }

    // RENEWAL "MOBILE": every component, the transported ones included, relaxes towards the reference solution
    // after the transported values are written (boundary compartments: the water of precipitation)
    if (OpenWQ_hostModelconfig.get_time_step() > 0 && !P.renewal_ref_c.empty() && !P.renewal_mobile.empty()) {
        const double dt = OpenWQ_hostModelconfig.get_time_step();
        for (int i = 0; i < nxyz; i++) {
            if (!active_cell(i) || P.renewal_ref[i] < 0 || !P.renewal_mobile[i]) continue;
            const double f = 1.0 - std::exp(-dt / P.renewal_tau_s[i]);
            const std::vector<double>& cref = P.renewal_ref_c[P.renewal_ref[i]];
            for (int j = 0; j < nc; j++) {
                if (P.is_diagnostic.size() > (size_t)j && P.is_diagnostic[j]) continue;
                c[(size_t)j * nxyz + i] += f * (cref[j] - c[(size_t)j * nxyz + i]);
            }
        }
    }

    // RESET_WATER (phreeqc_prepare_coupling): the H2O component back to its value at the first call
    if (P.reset_water) {
        if (P.h2o_ref.empty()) {
            P.h2o_ref.assign(nxyz, 0.0);
            for (int i = 0; i < nxyz; i++) P.h2o_ref[i] = c[(size_t)P.comp_H2O * nxyz + i];
        }
        static const bool dbg_h2o = std::getenv("OPENWQ_DEBUG_H2O") != nullptr;
        const long step = OpenWQ_hostModelconfig.get_interaction_step();
        for (int i = 0; i < nxyz; i++) {
            if (!active_cell(i)) continue;
            double& w = c[(size_t)P.comp_H2O * nxyz + i];
            if (dbg_h2o && i < 3 && (step <= 3 || step % 2000 == 0))
                std::fprintf(stderr, "H2O cell %d step %ld before %.10e ref %.10e\n", i, step, w, P.h2o_ref[i]);
            w = P.h2o_ref[i];
        }
    }

    indx = 0;

    // Validate and fix temperatures - PHREEQC needs reasonable values (0-100°C typical)
    for (int i = 0; i < nxyz; i++) {
        if (std::isnan(temperatures[i]) || std::isinf(temperatures[i]) || temperatures[i] > 200) {
            temperatures[i] = 25.0;  // Default to 25°C
        }
        if (std::isnan(pressures[i]) || std::isinf(pressures[i]) || pressures[i] <= 0) {
            pressures[i] = 1.0;  // Default to 1 atm
        }
        // Liquid water below 0 C (frozen soil layers, supercooled water) is
        // outside the range of the PHREEQC databases: the log K and gas
        // solubility expressions are fitted from 0 C up and extrapolate below
        // it (O2 solubility ~20 mg/L at -10 C). By default chemistry is
        // evaluated at 0 C instead. A compartment can lower this bound with
        // "MIN_C" in its TEMPERATURE entry, e.g. for the unfrozen water of
        // frozen soil, so that rate laws with a temperature factor see the
        // host temperature (rate laws that use a gas solubility should then
        // bound the temperature of that solubility themselves).
        // TEMPERATURE "SMOOTH_HOURS": exponential running mean of the mapped temperature (damped diel cycle)
        if (P.temp_tau_s[i] > 0.0 && OpenWQ_hostModelconfig.get_time_step() > 0) {
            double& m = P.temp_mean[i];
            if (!std::isfinite(m)) m = temperatures[i];
            else m += (1.0 - std::exp(-OpenWQ_hostModelconfig.get_time_step() / P.temp_tau_s[i])) * (temperatures[i] - m);
            {
                static const bool dbg = std::getenv("OPENWQ_DEBUG_TSM") != nullptr;
                const long step = OpenWQ_hostModelconfig.get_interaction_step();
                if (dbg && (step <= 24 || step % 2000 == 0))
                    std::fprintf(stderr, "TSM phreeqc %ld %.10e %.10e\n", step, temperatures[i], m);
            }
            temperatures[i] = m;
        }
        if (temperatures[i] < temperatures_min[i]) temperatures[i] = temperatures_min[i];
    }

    // Log first cell values for debugging (with units)
    if (num_chem_to_use > 0 && has_gfw) {
        double first_gfw = (gfw.size() > 0) ? gfw[0] : 1.0;
        msg_string = "<OpenWQ> PHREEQC: Setting T[0]=" + std::to_string(temperatures[0]) + " C"
            + ", P[0]=" + std::to_string(pressures[0]) + " atm"
            + ", c[0]=" + std::to_string(c[0]) + " mol/kgw"
            + " (GFW=" + std::to_string(first_gfw) + " g/mol)"
            + ", vol[0]=" + std::to_string(volumes[0]) + " L";
    } else {
        msg_string = "<OpenWQ> PHREEQC: Setting T[0]=" + std::to_string(temperatures[0])
            + ", P[0]=" + std::to_string(pressures[0])
            + ", c[0]=" + std::to_string(c[0])
            + ", vol[0]=" + std::to_string(volumes[0]);
    }
    if (OpenWQ_wqconfig.debug_mode) OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetPressure(pressures);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetTemperature(temperatures);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetConcentrations(c);
    // NOTE: Do NOT call SetRepresentativeVolume here!
    // The representative volume is set to 1.0 L during setup (phreeqc_setup.cpp).
    // We pass concentrations (mol/kgw) to PHREEQC, computed from mass/volume.
    // After PHREEQC returns updated concentrations, we convert back to mass
    // using the actual hydrological volumes. Calling SetRepresentativeVolume
    // with actual volumes (millions of liters) causes PHREEQC convergence failures.
    // OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetRepresentativeVolume(volumes);

    // Set fixed density and saturation to avoid calc_dens() crash
    std::vector<double> density(nxyz, 1.0);  // 1.0 kg/L
    std::vector<double> saturation(nxyz, 1.0);  // Fully saturated
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetDensityUser(density);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetSaturationUser(saturation);

    // CRITICAL: Tell PhreeqcRM to use user-provided density, not calculate internally
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->UseSolutionDensityVolume(false);

    msg_string = "<OpenWQ> PHREEQC: About to call RunCells with timestep="
        + std::to_string(OpenWQ_hostModelconfig.get_time_step()) + " s";
    if (OpenWQ_wqconfig.debug_mode) OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetTimeStep(OpenWQ_hostModelconfig.get_time_step());
    // Capture the RunCells return code (IRM_OK == 0). Previously the return was
    // ignored and "RunCells succeeded" was logged unconditionally, hiding
    // non-convergence failures (which then propagated as garbage/zeros).
    int phreeqc_runcells_status = 0;
    if (OpenWQ_hostModelconfig.get_time_step() != 0) {
        phreeqc_runcells_status = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->RunCells();
    }
    if (phreeqc_runcells_status != 0) {
        msg_string = "<OpenWQ> ERROR: PHREEQC RunCells did NOT converge (IRM_RESULT="
            + std::to_string(phreeqc_runcells_status)
            + "). Geochemistry this step is unreliable - check the .pqi solution, the "
            "mobile-species mapping and unit conversion (GFW).";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    } else {
        msg_string = "<OpenWQ> PHREEQC: RunCells succeeded";
        if (OpenWQ_wqconfig.debug_mode) OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

    // Get updated concentrations from PHREEQC (in mol/kgw)
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetConcentrations(c);
    msg_string = "<OpenWQ> PHREEQC: GetConcentrations succeeded, c.size()=" + std::to_string(c.size());
    if (OpenWQ_wqconfig.debug_mode) OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    // Read-out components from the speciated solutions (DIAGNOSTICS): pH as
    // -(log10 m + log10 gamma) of H+, divided by 1000 GFW so that the mg/L of
    // the component equal the pH; a species as its molality. Species absent
    // from a solution (concentration 0) give 0, not 10^0 (PhreeqcRM fills their
    // log10 molality with 0).
    if (!P.diagnostics.empty()) {
        std::vector<double> lm, lg, sc;
        P.phreeqcrm->GetSpeciesLog10Molalities(lm);
        P.phreeqcrm->GetSpeciesLog10Gammas(lg);
        P.phreeqcrm->GetSpeciesConcentrations(sc);
        const size_t nsp = P.phreeqcrm->GetSpeciesCount();
        const bool ok = lm.size() == nsp * (size_t)nxyz && lg.size() == lm.size() && sc.size() == lm.size();
        for (const auto& d : P.diagnostics) {
            double gfw_d = (has_gfw && (size_t)d.comp < gfw.size()) ? gfw[d.comp] : 1.0;
            if (!(gfw_d > 0.0) || !std::isfinite(gfw_d)) gfw_d = 1.0;
            for (int i = 0; i < nxyz; i++) {
                if (!active_cell(i)) continue;
                double v = 0.0;
                if (ok && d.species < 0) {
                    const size_t k = (size_t)P.species_Hplus * nxyz + i;
                    v = -(lm[k] + lg[k]) / (1000.0 * gfw_d);
                } else if (ok) {
                    const size_t k = (size_t)d.species * nxyz + i;
                    v = (sc[k] > 0.0) ? std::pow(10.0, lm[k]) : 0.0;
                }
                c[(size_t)d.comp * nxyz + i] = (std::isfinite(v) && v > 0.0) ? v : 0.0;
            }
        }
    }

    msg_string = "<OpenWQ> PHREEQC DEBUG: mobile_species.size()=" + std::to_string(OpenWQ_wqconfig.CH_model->PHREEQC->mobile_species.size())
        + ", num_chem=" + std::to_string(OpenWQ_wqconfig.CH_model->PHREEQC->num_chem)
        + ", num_chem_to_use=" + std::to_string(num_chem_to_use);
    if (OpenWQ_wqconfig.debug_mode) OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    // ########################################################
    // UNIT CONVERSION: PHREEQC (mol/kgw) -> OpenWQ (grams)
    // ########################################################
    // PHREEQC returns concentration in mol/kgw
    // Convert to g/L: conc_g_L = conc_mol_kgw * GFW
    // Convert to mass: mass_g = conc_g_L * volume_L
    // Delta mass = new_mass - old_mass
    // ########################################################

    for (unsigned int chemi=0;chemi<num_chem_to_use;chemi++){
        indx = 0;

        // Get the PHREEQC component index for this OpenWQ chemical
        unsigned int phreeqc_idx = OpenWQ_wqconfig.CH_model->PHREEQC->mobile_species[chemi];
        msg_string = "<OpenWQ> PHREEQC DEBUG: chemi=" + std::to_string(chemi)
            + ", phreeqc_idx=" + std::to_string(phreeqc_idx);
        if (OpenWQ_wqconfig.debug_mode) OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

        // Get GFW for this PHREEQC component
        double component_gfw = (has_gfw && phreeqc_idx < gfw.size()) ? gfw[phreeqc_idx] : 1.0;

        // Protect against invalid GFW
        if (component_gfw <= 0 || std::isnan(component_gfw) || std::isinf(component_gfw)) {
            component_gfw = 1.0;
        }

        for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){

            nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
            ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
            nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements
            for (int ix=0; ix<nx; ix++) {
                for (int iy = 0; iy<ny;iy++) {
                    for (int iz=0; iz<nz; iz++) {
                        // PHREEQC concentration in mol/kgw - use PHREEQC component index
                        double conc_mol_kgw = c[phreeqc_idx*nxyz+indx];

                        // Validate concentration
                        if (std::isnan(conc_mol_kgw) || std::isinf(conc_mol_kgw)) {
                            conc_mol_kgw = 0.0;
                        }
                        if (conc_mol_kgw < 0) {
                            conc_mol_kgw = 0.0;  // PHREEQC shouldn't return negative, but protect anyway
                        }

                        // Convert mol/kgw to g/L
                        // g/L = mol/kgw * GFW(g/mol)
                        double conc_g_L = conc_mol_kgw * component_gfw;

                        // Convert to mass: mass_g = conc_g_L * volume_L
                        double new_mass_g = conc_g_L * volumes[indx];

                        // Current mass in OpenWQ (index by PHREEQC component index)
                        double old_mass_g = (*OpenWQ_vars.chemass)(icmp)(phreeqc_idx)(ix,iy,iz);

                        // Delta mass for the chemistry derivative (same slot as transport).
                        // Cells that do not run PHREEQC (no reaction block) keep their mass.
                        const auto& cell_active = OpenWQ_wqconfig.CH_model->PHREEQC->cell_active;
                        const bool active = cell_active.empty() || cell_active[indx];
                        (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(phreeqc_idx)(ix,iy,iz) =
                            active ? (new_mass_g - old_mass_g) : 0.0;
                        if (OpenWQ_wqconfig.debug_mode && active) {
                            msg_string = "<OpenWQ> PHREEQC DEBUG cell " + std::to_string(indx)
                                + " comp " + std::to_string(phreeqc_idx)
                                + ": vol_L=" + std::to_string(volumes[indx])
                                + " old_g=" + std::to_string(old_mass_g)
                                + " c_out=" + std::to_string(conc_mol_kgw)
                                + " new_g=" + std::to_string(new_mass_g);
                            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                        }

                        indx++;
                    }
                }
            }
        }
    }

    // Log summary of first component change for debugging
    if (num_chem_to_use > 0 && nxyz > 0) {
        double first_gfw = (has_gfw && gfw.size() > 0) ? gfw[0] : 1.0;
        double conc_after_mol = c[0];
        double conc_after_g_L = conc_after_mol * first_gfw;
        msg_string = "<OpenWQ> PHREEQC: After reaction c[0]=" + std::to_string(conc_after_mol) + " mol/kgw"
            + " = " + std::to_string(conc_after_g_L) + " g/L"
            + " = " + std::to_string(conc_after_g_L * 1000) + " mg/L";
        if (OpenWQ_wqconfig.debug_mode) OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

}
