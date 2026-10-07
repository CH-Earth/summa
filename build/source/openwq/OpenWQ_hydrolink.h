// Copyright 2020, Diogo Costa (diogo.pinhodacosta@canada.ca)
// This file is part of OpenWQ model.

// This program, openWQ, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// =============================================================================
// OpenWQ hydrolink for SUMMA with (optional) internally coupled mizuRoute
// =============================================================================
// One OpenWQ instance holds the land compartments of SUMMA and, when the river
// network is coupled, the reaches of mizuRoute.
//
// Land compartments have one cell column per (HRU, domain) pair of SUMMA:
//   0 SCALARCANOPYWAT        canopy storage                    (1 layer)
//   1 ILAYERVOLFRACWAT_SNOW  snow layers                       (max snow layers)
//   2 RUNOFF                 water at the surface in the step  (1 layer, through-flow pool)
//   3 ILAYERVOLFRACWAT_SOIL  soil layers                       (max soil layers)
//   4 SCALARAQUIFER          aquifer                           (1 layer)
//   5 RUNOFF_TO_STREAM       water on its way to the stream    (1 layer, routing pool)
//   6 ILAYERVOLFRACWAT_LAKE  lake layers (only if any domain has lake layers)
// River compartment (only with mizuRoute), one cell per reach:
//   next index: RIVER_NETWORK_REACHES
// =============================================================================

#ifndef OPENWQ_HYDROLINK_INCLUDED
#define OPENWQ_HYDROLINK_INCLUDED

#include "global/OpenWQ_hostModelConfig.hpp"
#include "global/OpenWQ_json.hpp"
#include "global/OpenWQ_wqconfig.hpp"
#include "global/OpenWQ_vars.hpp"
#include "couplercalls/headerfile_CC.hpp"
#include "readjson/headerfile_nlohmann.hpp"
#include "initiate/headerfile_INIT.hpp"
#include "extwatflux_ss/headerfile_EWF_SS.hpp"
#include "units/headerfile_units.hpp"
#include "utils/headerfile_UTILS.hpp"
#include "compute/headerfile_compute.hpp"
#include "output/headerfile_OUT.hpp"
#include "models_CH/headerfile_CH.hpp"
#include "models_TD/headerfile_TD.hpp"
#include "models_LE/headerfile_LE.hpp"
#include "models_SI/headerfile_SI.hpp"
#include "models_TS/headerfile_TS.hpp"

#include <iostream>
#include <time.h>
#include <vector>
#include <filesystem>
#include <memory>
#include <string>

// Fixed compartment indices (the Fortran coupler uses the same numbers)
inline int canopy_index_openwq  = 0;
inline int snow_index_openwq    = 1;
inline int runoff_index_openwq  = 2;
inline int soil_index_openwq    = 3;
inline int aquifer_index_openwq = 4;
inline int stream_index_openwq  = 5;

// Flux-concentration exports (0-based, the Fortran coupler uses the same numbers)
inline int scalarRunoffVol_fluxexp_openwq     = 0;  // surface runoff [m3]
inline int averageRoutedRunoff_fluxexp_openwq = 1;  // routed runoff of the GRU [m3]
inline int scalarTotalRunoff_fluxexp_openwq   = 2;  // total runoff of the column [m3]
inline int reachOutflow_fluxexp_openwq        = 3;  // reach outflow [m3] (only with mizuRoute)

// Dependency variables available to the kinetic expressions
inline int sm_depend_openwq    = 0;  // SM          volumetric liquid water content [-]
inline int tair_depend_openwq  = 1;  // Tair_K      air temperature [K]
inline int tsoil_depend_openwq = 2;  // Tsoil_K     soil temperature [K]
inline int swrad_depend_openwq = 3;  // SWrad_Wm2   incoming shortwave radiation [W/m2]
inline int area_depend_openwq  = 4;  // cellArea_m2 area of the column or contributing area of the reach [m2]
inline int tc_depend_openwq    = 5;  // T           soil temperature (land) or air temperature (reach) [degC]

class CLASSWQ_openwq {

private:

    std::unique_ptr<OpenWQ_hostModelconfig> OpenWQ_hostModelconfig_ref =
        std::make_unique<OpenWQ_hostModelconfig>();
    std::unique_ptr<OpenWQ_couplercalls> OpenWQ_couplercalls_ref =
        std::make_unique<OpenWQ_couplercalls>();
    std::unique_ptr<OpenWQ_json> OpenWQ_json_ref =
        std::make_unique<OpenWQ_json>();
    std::unique_ptr<OpenWQ_wqconfig> OpenWQ_wqconfig_ref =
        std::make_unique<OpenWQ_wqconfig>();
    std::unique_ptr<OpenWQ_units> OpenWQ_units_ref =
        std::make_unique<OpenWQ_units>();
    std::unique_ptr<OpenWQ_utils> OpenWQ_utils_ref =
        std::make_unique<OpenWQ_utils>();
    std::unique_ptr<OpenWQ_readjson> OpenWQ_readjson_ref =
        std::make_unique<OpenWQ_readjson>();
    std::unique_ptr<OpenWQ_initiate> OpenWQ_initiate_ref =
        std::make_unique<OpenWQ_initiate>();
    std::unique_ptr<OpenWQ_extwatflux_ss> OpenWQ_extwatflux_ss_ref =
        std::make_unique<OpenWQ_extwatflux_ss>();
    std::unique_ptr<OpenWQ_compute> OpenWQ_solver_ref =
        std::make_unique<OpenWQ_compute>();
    std::unique_ptr<OpenWQ_output> OpenWQ_output_ref =
        std::make_unique<OpenWQ_output>();
    std::unique_ptr<OpenWQ_vars> OpenWQ_vars_ref;
    std::unique_ptr<OpenWQ_TD_model> OpenWQ_TD_model_ref =
        std::make_unique<OpenWQ_TD_model>();
    std::unique_ptr<OpenWQ_LE_model> OpenWQ_LE_model_ref =
        std::make_unique<OpenWQ_LE_model>();
    std::unique_ptr<OpenWQ_CH_model> OpenWQ_CH_model_ref =
        std::make_unique<OpenWQ_CH_model>();
    std::unique_ptr<OpenWQ_SI_model> OpenWQ_SI_model_ref =
        std::make_unique<OpenWQ_SI_model>();
    std::unique_ptr<OpenWQ_TS_model> OpenWQ_TS_model_ref =
        std::make_unique<OpenWQ_TS_model>();

    int num_col   = 0;            // land columns: (HRU, domain) pairs
    int num_reach = 0;            // river reaches (0 without mizuRoute)
    int nz_depend = 1;            // layers of the dependency variables
    int lake_index_openwq  = -1;  // lake compartment (-1 if absent)
    int river_index_openwq = -1;  // river compartment (-1 if absent)

    time_t to_time(int simtime_summa[]);

public:

    CLASSWQ_openwq();
    ~CLASSWQ_openwq();

    // Declare compartments, external fluxes, exports and dependencies, and read the configuration
    int decl(
        int num_col,
        int nCanopy_2openwq,
        int nSnow_2openwq,
        int nSoil_2openwq,
        int nRunoff_2openwq,
        int nAquifer_2openwq,
        int nLake_2openwq,
        int nYdirec_2openwq,
        long long colId[],
        int colDom[],
        int has_glacier,
        int num_reach,
        long long reachId[]);

    // Water volume and dependencies of the reaches at the start of the step (call before the columns)
    int openwq_set_reach_state(
        int n_reach,
        double reachVol_m3[],
        double airTemp_K[],
        double SWrad_Wm2[],
        double area_m2[]);

    // Water volumes and dependencies of one land column at the start of the step;
    // the last column starts the OpenWQ time step
    int openwq_run_time_start(
        int index_col,
        int nSnow_2openwq,
        int nLake_2openwq,
        int nSoil_2openwq,
        int simtime_summa[],
        double soilMoist_depVar_summa_frac[],
        double soilTemp_depVar_summa_K[],
        double airTemp_depVar_summa_K,
        double SWrad_depVar_summa_Wm2,
        double sweWatVol_stateVar_summa_m3[],
        double lakeWatVol_stateVar_summa_m3[],
        double canopyWatVol_stateVar_summa_m3,
        double soilWatVol_stateVar_summa_m3[],
        double aquiferWatVol_stateVar_summa_m3,
        double col_area_m2);

    // Water flux between two cells (indices are 1-based; recipient -1 leaves the domain).
    // out_to_stream = 1: the solute that leaves the domain is collected in RUNOFF_TO_STREAM.
    int openwq_run_space(
        int simtime_summa[],
        int source, int ix_s, int iy_s, int iz_s,
        int recipient, int ix_r, int iy_r, int iz_r,
        double wflux_s2r, double wmass_source,
        int out_to_stream);

    // Water flux entering the domain from an external source (indices are 1-based)
    int openwq_run_space_in(
        int simtime_summa[],
        std::string source_EWF_name,
        int recipient, int ix_r, int iy_r, int iz_r,
        double wflux_s2r);

    // Water volume of one cell (indices are 1-based)
    int openwq_set_watervol(
        int icmp, int ix, int iy, int iz, double vol_m3);

    // Through-volume of a flux-concentration export (iflux is 0-based, cell indices are 1-based)
    int openwq_set_fluxvol(
        int iflux, int ix, int iy, int iz, double flux_vol_m3);

    // Solve the time step and write the output
    int openwq_run_time_end(
        int simtime_summa[]);
};

#endif // OPENWQ_HYDROLINK_INCLUDED
