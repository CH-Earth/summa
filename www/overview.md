**SUMMA** or the **Structure for Unifying Multiple Modeling Alternatives** is a hydrologic modeling aproach that is built on a common set of governing equations and a common numerical solver, which together constitute the “structural core” of the model. Different modeling approaches can then be implemented within the structural core, enabling a controlled and systematic analysis of alternative modeling options, and providing insight for future model development.

The important modeling features are:

 1. The formulation of the governing model equations is cleanly separated from their numerical solution;

 1. Different model representations of physical processes (in particular, different flux parameterizations) can be used within a common set of governing equations; and

 1. The physical processes can be organized in different spatial configurations, including model elements of different shape and connectivity (e.g., nested multi-scale grids and HRUs).

SUMMA can be used to configure a wide range of hydrological model alternatives. We anticipate that systematic model analysis will help researchers and practitioners understand reasons for inter-model differences in model behavior, and, when applied across a large sample of catchments, may provide insights on the dominance of different physical processes and regional variability in the suitability of different modeling approaches. An important application of SUMMA is selecting specific physics options to reproduce the behavior of existing models – these applications of *“model mimicry”* can be used to define reference (benchmark) cases in structured model comparison experiments, and can help diagnose weaknesses of individual models in different hydroclimatic regimes.