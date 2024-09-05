// This program, openWQ, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef OPENWQ_HYDROLINK_INCLUDED
#define OPENWQ_HYDROLINK_INCLUDED

#include "couplercalls/OpenWQ_couplercalls.hpp"
#include "global/OpenWQ_hostModelconfig.hpp"
#include "global/OpenWQ_json.hpp"
#include "global/OpenWQ_wqconfig.hpp"
#include "global/OpenWQ_vars.hpp"
#include "readjson/OpenWQ_readjson.hpp"
#include "initiate/OpenWQ_initiate.hpp"
#include "chem/OpenWQ_chem.hpp"
#include "watertransp/OpenWQ_watertransp.hpp"
#include "extwatflux_ss/OpenWQ_extwatflux_ss.hpp"
#include "units/OpenWQ_units.hpp"
#include "utils/OpenWQ_utils.hpp"
#include "solver/OpenWQ_solver.hpp"
#include "output/OpenWQ_output.hpp"
#include <iostream>
#include <time.h>
#include <vector>
#include <filesystem>

// Global Indexes for Compartments
  inline int canopy_index_openwq    = 0;
  inline int snow_index_openwq      = 1;
  inline int runoff_index_openwq    = 2;
  inline int soil_index_openwq      = 3;
  inline int aquifer_index_openwq   = 4;
  inline int max_snow_layers        = 5;

class CLASSWQ_openwq
{

    // Instance Variables
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
        std::unique_ptr<OpenWQ_watertransp> OpenWQ_watertransp_ref = 
            std::make_unique<OpenWQ_watertransp>();
        std::unique_ptr<OpenWQ_chem> OpenWQ_chem_ref = 
            std::make_unique<OpenWQ_chem>();
        std::unique_ptr<OpenWQ_extwatflux_ss> OpenWQ_extwatflux_ss_ref = 
            std::make_unique<OpenWQ_extwatflux_ss>();
        std::unique_ptr<OpenWQ_solver> OpenWQ_solver_ref = 
            std::make_unique<OpenWQ_solver>();
        std::unique_ptr<OpenWQ_output> OpenWQ_output_ref = 
            std::make_unique<OpenWQ_output>();

        std::unique_ptr<OpenWQ_vars> OpenWQ_vars_ref; // Requires input from summa 

        int num_HRU;
        const float *hru_area;

    // Constructor
    public:
        CLASSWQ_openwq();
        ~CLASSWQ_openwq();
    
    // Methods
    void printNum() {
        std::cout << "num = " << this->num_HRU << std::endl;
    }

    int decl(
        int num_HRU,                // num HRU
        int nCanopy_2openwq,      // num layers of canopy (fixed to 1)
        int nSnow_2openwq,        // num layers of snow (fixed to max of 5 because it varies)
        int nSoil_2openwq,        // num layers of snoil (variable)
        int nRunoff_2openwq,      // num layers of runoff (fixed to 1)
        int nAquifer_2openwq,     // num layers of aquifer (fixed to 1)
        int nYdirec_2openwq);           // num of layers in y-dir (set to 1 because not used in summa)

    int openwq_run_time_start(
        bool last_hru_flag,
        int hru_index, 
        int nSnow_2openwq, 
        int nSoil_2openwq, 
        int simtime_summa[],
        double soilMoist_depVar_summa_frac[],                  
        double soilTemp_depVar_summa_K[],
        double airTemp_depVar_summa_K,
        double sweWatVol_stateVar_summa_m3[],
        double canopyWatVol_stateVar_summa_m3,
        double soilWatVol_stateVar_summa_m3[],
        double aquiferWatVol_stateVar_summa_m3);

    int openwq_run_space(
        int simtime_summa[], 
        int source, int ix_s, int iy_s, int iz_s,
        int recipient, int ix_r, int iy_r, int iz_r, 
        double wflux_s2r, double wmass_source);

    int openwq_run_space_in(
        int simtime_summa[],
        std::string source_EWF_name,
        int recipient, int ix_r, int iy_r, int iz_r, 
        double wflux_s2r);

    int openwq_run_time_end(
        int simtime_summa[]);
        
    int get_numHRU();
        
};
#endif