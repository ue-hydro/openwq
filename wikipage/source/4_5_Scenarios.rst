Scenario Analysis
==================================

OpenWQ has a dedicated scenario folder, ``supporting_scripts/4_Scenarios``,
built like the calibration one: the **scenario template**
(``scenario_config_template.py``) and the ``scenario_lib`` package (the option
catalogue and runner ``scenarios.py``, the setup report
``Gen_Scenario_Setup_Report.py`` and the comparison report
``Gen_Scenario_Results_Report.py``), which builds on
``3_Calibration/calibration_lib``. It serves *what-if* analysis on the
**calibrated model**: management practices, land-use change, point sources and wastewater
treatment, governance targets, climate change, in-stream measures and legacy
nutrient stores. Every scenario is a named list of **options** (levers) of two
groups — **Climate** and **Management & policy** — and every option is
expressed through the model's **input files** (source/sink loads, initial
conditions, a biogeochemical rate map, the host-model forcing) — nothing in the
engine changes — so the same catalogue works for every host model, load source
and module combination.

The scenario runner copies the calibrated model's inputs into one run folder per
scenario, applies the levers, runs the coupled model once per scenario plus the
*baseline* (the calibrated model with no lever), and writes a comparison report
in which **the baseline and every scenario appear on the same graphs**.


Workflow
~~~~~~~~~~~~~~~~~~~

The scenario template sits at the end of the usual chain:

1. **Model config** — run the model-config template, check its config report and
   use its snippets to run the model and produce the output report.
2. **Calibration** — run the calibration template, check the calibration config
   report, start the calibration from its snippets and, when it finishes, open
   the calibration results report.
3. **Save the calibrated model config** — the results report has a box *"Save
   the calibrated model config"* that writes ``<model_template>_config_run.py``
   next to your model-config template. It is your own template plus a short
   block that (a) writes the inputs to ``<dir2save_input_files>_calibrated``,
   (b) bakes the generation-time calibrated values (climate-response parameters
   of the loads) and (c) tells the input generator to apply the calibrated best
   parameter values and the ML layers to the freshly generated **full-period**
   inputs, exactly as the best evaluation was configured. Run it once like any
   model-config template::

       python <model_template>_config_run.py

   It produces the calibrated inputs (+ ``calibrated_setup.json`` with the
   provenance: calibration folder, best evaluation, objective, applied values)
   and the usual config report with the snippets to run the calibrated model
   and to build its output report.
4. **Scenario template** — copy ``4_Scenarios/scenario_config_template.py``
   next to the case, point ``config_run_path`` at the calibrated config and run
   it::

       python scenario_config_template.py

   The **scenario config report** opens (work dir
   ``<dir2save_input_files>_calibrated_scenarios`` by default):

   * *Overview* tab — the base model and its calibrated setup (a warning is shown
     when the base has no ``calibrated_setup.json``, i.e. it is not calibrated,
     or has not been run yet);
   * *Scenarios* tab — the **scenario set**, one tab per scenario. *+ New
     scenario* adds an empty scenario that you define yourself; add as many
     as needed. Each scenario is edited in two blocks:

     - **Climate** (whole domain): tick *Climate change* (delta-change) and /
       or *Precipitation intensification* and set the changes.
     - **Management & policy**, organised by **where it applies**:

       - **All units** (HRUs, sub-basins or reaches): the options ticked here
         apply everywhere.
       - **Groups of units**: *+ Group of …* creates a group; pick its units
         on the map (click to add / remove, *Done*) or type their ids, then
         tick the group's own options. A unit belongs to one group. In a
         group, an option **replaces** the same option set for all units (for
         example a stronger fertilizer reduction in the headwaters); different
         options add up. Any number of groups can be added, and each scenario
         has its own.

       Every zone offers one collapsible block per category — *Land-use
       change*, *Agricultural practices*, *Urban areas & point sources*,
       *Policy targets & regulation*, *In-stream measures & legacy stores*,
       *Generic adjustment*. Tick an option to open its inputs. Where a
       published relation exists (see *Estimating an option from
       measurements* below) the inputs are **measurable quantities** — strip
       width, field slope, application rates, soil-test P, treatment class,
       animal numbers — and the change of the N and P export is **computed
       and shown** next to the option, with the relation and its reference;
       the computed percentages can still be overridden under *details*.
       Options without such a relation show their documented default
       percentages for adjustment, with the source. Repeatable options (a
       land-use conversion, a new point source, an in-stream rate change, a
       load multiplier) are added with their *+ Add …* button.

     The map colours the groups and, on hover, lists the options that apply
     to a unit. It is drawn from the model config's basin / river-network
     shapefile; the id column is detected from the unit ids of the run.
     Without a readable shapefile the ids are typed.

     A **Summary** table lists, per scenario, the number of climate options,
     of options for all units and the groups, and whether anything is missing.
     Then set the **run settings**: scored period (defaults to the model
     period), spin-up start, parallel runs, species / compartment to compare,
     concentration thresholds for exceedance statistics, and the container
     runtime (Docker locally, Apptainer on HPC);
   * right pane — the generated scenario script (updated live) with the
     commands to run it and where the results land.
5. **Run the scenarios** — save the scenario script (right pane, step 1) and
   run it::

       python <template>_scenarios_run.py            # runs baseline + scenarios
       python <template>_scenarios_run.py --report   # rebuilds the report only

The scenario script lists the scenarios as plain Python data and can be edited
by hand::

    scenarios = [
        {"name": "headwater_measures",
         "climate": [ {"id": "climate_delta", "params": {"dP_pct": "4", "dT_c": "1.5"}} ],
         "management": [
             {"where": "all",
              "options": [ {"id": "fert_reduction", "params": {"reduction_N_pct": 30}} ]},
             {"name": "headwaters", "where": ["12", "15", "16"],
              "options": [ {"id": "fert_reduction", "params": {"reduction_N_pct": 50}},
                           {"id": "riparian_buffer", "params": {"adoption_pct": 100}} ]},
         ]},
    ]

Outputs go to the scenario work dir: ``scenarios/<name>/`` — one run folder per
scenario (a complete, re-runnable openWQ setup with ``scenario_applied.json``
listing every lever and what it changed), ``scenarios_simulated.csv``,
``scenarios_summary.csv`` and the comparison report
``<template>_scenarios_report.html``.

.. note::

   Any model-config template can be the base of a scenario analysis (the
   scenario template only needs a config whose inputs exist), but the point is
   the calibrated one: the ``_config_run.py`` saved from the results report.


The option catalogue
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The catalogue lives in ``4_Scenarios/scenario_lib/scenarios.py`` (``GROUPS``,
``CATEGORIES``, ``LEVERS``) and is the single source of truth for the tab and
the runner: 50 options. It was checked against the management and scenario
options of SWAT / SWAT+, HYPE, INCA, the MIKE family (MIKE HYDRO Basin and its
Load Calculator, MIKE SHE-DAISY), MONERIS, GWLF-E / MapShed, the Chesapeake Bay
Program BMP list, HSPF, APEX, AnnAGNPS and SPARROW. Options on the diffuse
loads use the per-class load breakdown the load generator writes to
``openwq_in/ss_copernicus_files/nutrient_loads_detailed.csv`` (unit × land-use
class × species), so a class-based practice only touches the share of each
unit's load that comes from those classes, in the units of its zone. Options
compose multiplicatively.

**Climate** (whole domain)

* *Climate change* — delta-change on the host forcing (precipitation factor,
  temperature shift; annual or 12 monthly values) with the load generator's
  climate response ``precip^power × Q10^(ΔT/10)``; under *details*, changes of
  shortwave radiation, humidity and wind speed (SUMMA forcing). A perturbed
  copy of the forcing (SUMMA forcing NetCDF; mizuRoute runoff input scaled as
  a proxy) is written under ``forcing_scenario/`` in the run folder and the
  per-run file manager / control points at it. Presets are illustrative
  global-mean ranges — use regional CMIP6 deltas. A CO2 change is not
  represented. With SUMMA internally coupled to mizuRoute the single run
  carries the land and the river compartments (the executable receives the
  mizuRoute TOML with ``-c``); the compartment compared by default is then
  ``RIVER_NETWORK_REACHES`` and the statistics are summarised per compartment.
* *Precipitation intensification* — the wettest days are scaled up and the
  rest rescaled so the annual total is unchanged.

**Management & policy** (for all units, or for a group of units)

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Category
     - Options
   * - Land-use change & disturbance
     - Land-use conversion; Forest harvesting; Wildfire / prescribed burning. A conversion recomputes the unit's load from the
       per-class export coefficients (presets: afforestation, reforestation
       of pasture, land retirement, perennial energy crops, agricultural
       expansion, urbanization); harvesting and fire *increase* the export of
       the affected share of the classes.
   * - Agricultural practices
     - 29 practices (table below).
   * - Urban areas & point sources
     - Urban stormwater BMPs; Urban nutrient management; Construction-site erosion and sediment control; Septic systems: upgrade or connection to sewers; Sewer system: overflow storage and rehabilitation; Scale existing point-source loads; New point source. Urban BMP presets follow the Chesapeake Bay
       Program efficiencies (wet ponds, detention, infiltration, filtering,
       bioretention, open channels, permeable pavement, street sweeping). A
       new point source is a constant kg/day at one unit, optionally growing
       every year.
   * - Policy targets & regulation
     - Load-reduction target; Maximum application rate per hectare; Wastewater treatment standard; Atmospheric nitrogen deposition change; Phosphate-free detergents.
   * - In-stream measures & legacy stores
     - Stream restoration / floodplain reconnection; Enhanced in-stream processing; Legacy nutrient stores. In-stream processing multiplies one biogeochemical
       rate constant in the zone's units only (the Layer-1 spatial-parameter
       machinery, see :doc:`4_4_Hybrid_ML`); legacy stores scale the initial
       conditions (whole domain).
   * - Generic adjustment
     - Free-form load multiplier by species, class, unit and month.

**Agricultural practices.** Each practice cuts the nutrient export of the
selected land-use classes by a percentage for the N species, the P species and
(under *details*) every other species, on a share of the class area
(*adoption*), in the units of its zone, optionally from a *start year* on
(phased adoption). A negative value increases the export (dissolved P under
no-till, residue removal, new tile drains). Defaults are mid-range values from
the Iowa Nutrient Reduction Strategy science assessment, the Chesapeake Bay
Program BMP efficiencies (2011 table), the SWAT documentation and the Swedish
programme of measures; where no documented range exists the option's
description says the default is illustrative. Replace them with regional
numbers.

.. list-table::
   :header-rows: 1
   :widths: 56 14 30

   * - Practice
     - Classes
     - Default export reduction
   * - **Nutrient & manure management**
     -
     -
   * - Fertilizer / manure rate reduction
     - crop
     - N 30 %, P 30 %
   * - Soil-test-based phosphorus application
     - crop
     - N 0 %, P 17 %
   * - Application timing / spreading window
     - crop
     - closed months
   * - Fertilizer placement & enhanced-efficiency products
     - crop
     - N 10 %, P 25 %
   * - Manure management
     - crop
     - N 20 %, P 20 %
   * - Animal feed management
     - grass
     - N 5 %, P 15 %
   * - **Crop & soil management**
     -
     -
   * - Cover crops / catch crops
     - crop
     - N 31 %, P 10 %
   * - Crop rotation / diversification
     - crop
     - N 20 %, P 10 %
   * - Conservation tillage / no-till
     - crop
     - N 10 %, P 30 %
   * - Tillage timing
     - crop
     - N 10 %, P 0 %
   * - Residue management
     - crop
     - N 5 %, P 20 %
   * - Irrigation management
     - crop
     - N 15 %, P 10 %
   * - Organic / low-input farming
     - crop
     - N 20 %, P 10 %
   * - **Erosion & runoff control**
     -
     -
   * - Contour farming
     - crop
     - N 10 %, P 30 %
   * - Strip cropping
     - crop
     - N 20 %, P 35 %
   * - Terracing
     - crop
     - N 20 %, P 50 %
   * - Grassed waterways
     - crop
     - N 15 %, P 35 %
   * - Water and sediment control basins / check dams
     - crop
     - N 0 %, P 85 %
   * - **Edge of field & drainage**
     -
     -
   * - Vegetative filter strips
     - crop
     - N 30 %, P 40 %
   * - Riparian buffer zones
     - all
     - N 40 %, P 40 %
   * - Drainage water management
     - crop
     - N 33 %, P 0 %
   * - Constructed wetlands / retention ponds
     - all
     - N 40 %, P 40 %
   * - **Livestock**
     -
     -
   * - Grazing management / livestock exclusion
     - grass
     - N 10 %, P 25 %
   * - Livestock numbers
     - grass
     - N 20 %, P 20 %
   * - Barnyard / feedlot runoff control
     - grass
     - N 20 %, P 20 %, other 40 %
   * - **Soil amendments**
     -
     -
   * - Biochar soil amendment
     - crop
     - N 13 %, P 0 %
   * - Structure liming of clay soils
     - crop
     - N 0 %, P 30 %
   * - **Other**
     -
     -
   * - User-defined conservation practice
     - crop
     - N 0 %, P 0 %
   * - Pesticide / other contaminant application reduction
     - crop
     - N 0 %, P 0 %, other 30 %

Options that cannot act on a run (a treatment upgrade when the run has no
point-source entries, a practice on a class that does not exist in the basin,
a land-use conversion when the load breakdown is absent, a unit with no loads)
are reported as such in the Summary table and in the results report; they
never fail the run silently.

**Not covered.** Options of other models that change the *water balance* of
the host model cannot be expressed through the openWQ input files: reservoir
and dam operation, water abstraction and allocation, inter-basin transfers,
paddy / impoundment water management, and the CO2 effect on plants. Run the
host model with the changed hydrology and use it as the base of the scenarios.


Estimating an option from measurements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the options below the report computes the change of the N and P export
from quantities the user can measure or look up, with a published relation.
The relation and its reference are shown next to the result; the export
coefficients of the run (per land-use class and species, from the load
generator) are used where the relation needs them. "Read" means the equation
was checked in the source itself; "reported" means it was confirmed through
a secondary source only.

.. list-table::
   :header-rows: 1
   :widths: 22 24 40 14

   * - Option
     - Inputs
     - Relation and reference
     - Checked
   * - Vegetative filter strips
     - strip width (m)
     - removal = K·(1 − e\ :sup:`−b·w`); K, b = 92.0, 0.160 (N), 89.5, 0.157 (P),
       90.9, 0.446 (sediment). Zhang et al. 2010, J. Environ. Qual. 39:76–84
       (meta-analysis of 73 studies, R² 0.44 for N, 0.35 for P).
     - read
   * - Riparian buffer zones
     - width (m), vegetation
     - N: 50, 75, 90 % removal at 17, 51, 84 m (herbaceous) or 3, 18, 44 m
       (forest / mixed), interpolated (Mayer et al. 2007, J. Environ. Qual.
       36:1172–1180); P: the Zhang et al. 2010 curve.
     - reported / read
   * - Constructed wetlands, retention ponds
     - wetland area, drained area, annual runoff, areal rate constants
     - first-order areal model, removal = 1 − e\ :sup:`−k/q` with q = runoff ×
       drained area / wetland area (Kadlec & Wallace 2009, k-C* model); k for TN
       median 12.6 m/yr over 116 free-water-surface wetlands; k for P a typical
       value to adjust.
     - reported
   * - Drainage water management
     - measure, share of the drain flow reduced, diverted or treated, nitrate
       removal in that share
     - N change = share × removal. Controlled drainage: the nitrate load follows
       the drain flow (Ross et al. 2016; meta-analysis mean −50 %, 19–82 %);
       saturated buffers: mean −44 ± 26 % of the field load, driven by the
       diverted share (Jaynes & Isenhart 2019); bioreactors: removal of the
       treated flow set by retention time and temperature, Q10 2.15 (Addy et
       al. 2016).
     - reported
   * - Contour farming, strip cropping, terracing
     - field slope, rotation or outlet type, particulate share of the N and P
       export
     - USLE support-practice factor P by slope class (contouring 0.50–0.90,
       strip cropping 0.25–0.90, terraces 0.05–0.18); the sediment-bound share
       of the export falls by 1 − P. Wischmeier & Smith 1978, USDA Agriculture
       Handbook 537, tables 13–15.
     - read
   * - Fertilizer / manure rate
     - N and P applied now and after, leached fraction, share of the P export
       from applied P
     - N change = FracLEACH × rate change / export coefficient of the classes;
       FracLEACH = 0.24 (0.01–0.73) in wet climates, 0 in dry climates (IPCC
       2019 Refinement, Vol. 4 Ch. 11, table 11.3). P change proportional to
       the P rate change.
     - read
   * - Soil-test-based P application
     - soil-test P now and target, dissolved share of the P export
     - dissolved P in runoff proportional to soil-test P (extraction
       coefficients 1.2–3.0 for Mehlich-3 / Bray-1); the dissolved share
       changes by the ratio of soil-test P values. Vadas et al. 2005, J.
       Environ. Qual. 34:572–580.
     - read
   * - Wastewater treatment standard
     - treatment class now and after
     - removal by class: primary 10 % N and P; secondary 35 % N, 45 % P;
       tertiary 80 % N, 90 % P; load change = 1 − (1 − removal after) /
       (1 − removal now). Van Drecht et al. 2009, Global Biogeochem. Cycles 23.
     - reported
   * - Septic systems
     - people upgraded, share of failing systems, measure
     - 12.0 g N and 2.5 g P per person per day in the effluent, 1.6 g N and
       0.4 g P per day taken up by plants in the growing season, no P from a
       normal system (GWLF manual, Haith et al. 1992, table B-18); connection
       removes the load, denitrifying units 50 % N, pumping 5 % N (CBP 2014);
       relative to the urban-class export of the basin.
     - read / reported
   * - Phosphate-free detergents
     - wastewater share of the urban P export
     - per-person P falls from 2.5 to 1.5 g/day (−40 %) (GWLF manual, table
       B-18).
     - read
   * - Livestock numbers
     - animals now and after, livestock share of the class export
     - excretion scales with animal numbers (IPCC 2019, Vol. 4 Ch. 10).
     - reported
   * - Animal feed management
     - share of animals on phytase, cut of P excretion, manure share of the P
       export
     - phytase lowers P excretion by 15–30 %.
     - reported
   * - Atmospheric N deposition
     - deposition now and after, fraction reaching the river
     - N change = fraction × deposition change / export coefficient; net
       anthropogenic N input studies find 20–25 % of the inputs exported by
       rivers (Howarth et al. 1996, 2012).
     - reported

Other relations that exist but are not (yet) wired to an option: the SWAT
vegetative-filter-strip model (White & Arnold 2009; SWAT theoretical
documentation ch. 6:5: runoff reduction RR = 75.8 − 10.8 ln RL + 25.9 ln
K\ :sub:`sat`, TN = 0.036 SR\ :sup:`1.69`, nitrate = 39.4 + 0.584 RR, TP =
0.90 SR, soluble P = 29.3 + 0.51 RR), the Simple Method for urban loads
(Schueler 1987: L = 0.226 R C A, R\ :sub:`v` = 0.05 + 0.9 I\ :sub:`a`) and
the Chesapeake retrofit curves for urban BMPs, the RivR-N in-stream retention
relation (Seitzinger et al. 2002: N removed = 88.45 (depth / travel
time)\ :sup:`−0.3677`), the SWAT pond settling model (mass settled =
v·c·A·dt, apparent settling velocity), and the Iowa CREP wetland
removal-versus-loading curves.

**Documented mean only.** Cover crops (Thapa et al. 2018: −56 % nitrate
leaching for non-legumes, more with early sowing and biomass), nitrification
inhibitors, conservation tillage and residue management, crop rotation,
biochar, structure liming, grazing measures, barnyard runoff control, urban
nutrient management, construction erosion control, street sweeping, sewer
measures, stream restoration (the Chesapeake protocols need bank-erosion
surveys and reach lengths), wildfire (Smith et al. 2011: first-year nitrate
export 3–250 times the unburnt value, scaling with the burnt area of moderate
and high severity). These options show their default percentages with the
source.

**No usable relation found.** Tillage timing, organic or low-input farming,
forest harvesting, fertilizer timing and placement, pesticide reduction.


Reading the results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The comparison report puts the **baseline (the calibrated model, no lever) on
every graph next to the scenarios**: the mean concentration per species
(absolute, baseline included), the **% change of the mean concentration versus
the baseline** (basin aggregate and per unit), the **threshold exceedance**
(share of the scored time steps above a threshold), and **time-series
overlays** for the observation units first, then the units with the largest
baseline concentrations. The recipe table lists, per scenario and option, its
group, its parameters, where it was applied and what the runner actually
changed. Every scenario's run folder is a complete,
re-runnable openWQ setup.

.. note::

   Scenarios change the *inputs* of a calibrated model; the calibrated
   parameters themselves are kept fixed. Whether a measure's effect is
   realistic therefore depends on the model's own sensitivity — a rate
   constant that does nothing in the calibrated model will do nothing in a
   scenario either.
