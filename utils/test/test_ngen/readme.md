# SUMMA case studies
This folder contains a case study to show how a typical SUMMA setup looks in the NextGen setup.The folder serves a double purpose as a way to track default versions of certain input files, such as the Noah-MP tables and spatially constant parameter files. These files are in /settings/SUMMA.


## Settings
The SUMMA folder contains the setting some files that typically do not change for different model applications. Currently these include:
- `TBL_GENPARM.TBL`: lookup table for general parameters (legacy, currently unused)
- `TBL_MPTABLE.TBL`: lookup table for vegetation parameters
- `TBL_SOILPARM.TBL`: lookup table for soil parameters
- `TBL_VEGPARM.TBL`: lookup table for vegetation parameters

## Routing
`domain_provo` and `gauge_01073000` include realization configs that route through
ngen with t-route (`*_routing.json` / `provo_routing.yaml`), reading the real
hydrofabric network directly from the domain's `.gpkg` file. For routing the
same domains with SUMMA's built-in mizuRoute coupling instead (either a
synthetic network, or the same real hydrofabric network for comparison against
a t-route run), see [`../test_mizuroute/`](../test_mizuroute/README.md).
