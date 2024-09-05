// Copyright 2020, Diogo Costa, diogo.pinhodacosta@canada.ca
// This file is part of OpenWQ model.

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
#include "OpenWQ_hydrolink.h"
#include "OpenWQ_interface.h"

// Constructor
// initalize numHRUs value
CLASSWQ_openwq::CLASSWQ_openwq() {}

// Deconstructor
CLASSWQ_openwq::~CLASSWQ_openwq() {}

int CLASSWQ_openwq::decl(
    int num_HRU,              
    int nCanopy_2openwq,      // num layers of canopy (fixed to 1)
    int nSnow_2openwq,        // num layers of snow (fixed to max of 5 because it varies)
    int nSoil_2openwq,        // num layers of snoil (variable)
    int nRunoff_2openwq,      // num layers in the runoff of SUMMA
    int nAquifer_2openwq,     // num layers of aquifer (fixed to 1)
    int nYdirec_2openwq){     // num of layers in y-dir (set to 1 because not used in summa)

    this->num_HRU = num_HRU;

    if (OpenWQ_hostModelconfig_ref->get_num_HydroComp()==0) {

        // Compartment names
        // Make sure to use capital letters for compartment names
        OpenWQ_hostModelconfig_ref->add_HydroComp(canopy_index_openwq,"SCALARCANOPYWAT", num_HRU, nYdirec_2openwq, nCanopy_2openwq);      // Canopy
        OpenWQ_hostModelconfig_ref->add_HydroComp(snow_index_openwq,"ILAYERVOLFRACWAT_SNOW", num_HRU, nYdirec_2openwq, max_snow_layers);  // snow (layerd)
        OpenWQ_hostModelconfig_ref->add_HydroComp(runoff_index_openwq,"RUNOFF", num_HRU, nYdirec_2openwq, nRunoff_2openwq);               // Runoff
        OpenWQ_hostModelconfig_ref->add_HydroComp(soil_index_openwq,"ILAYERVOLFRACWAT_SOIL", num_HRU, nYdirec_2openwq, nSoil_2openwq);    // Soil (layerd)
        OpenWQ_hostModelconfig_ref->add_HydroComp(aquifer_index_openwq,"SCALARAQUIFER", num_HRU, nYdirec_2openwq, nAquifer_2openwq);      // Aquifer

        // External fluxes
        // Make sure to use capital letters for external fluxes
        OpenWQ_hostModelconfig_ref->add_HydroExtFlux(0,"PRECIP", num_HRU,nYdirec_2openwq,1);


        OpenWQ_vars_ref = std::make_unique<OpenWQ_vars>(
            OpenWQ_hostModelconfig_ref->get_num_HydroComp(),
            OpenWQ_hostModelconfig_ref->get_num_HydroExtFlux());

        
        // Dependencies
        // to expand BGC modelling options
        OpenWQ_hostModelconfig_ref->add_HydroDepend(0,"SM",        num_HRU,nYdirec_2openwq, nSnow_2openwq + nSoil_2openwq);
        OpenWQ_hostModelconfig_ref->add_HydroDepend(1,"Tair_K",    num_HRU,nYdirec_2openwq, nSnow_2openwq + nSoil_2openwq);
        OpenWQ_hostModelconfig_ref->add_HydroDepend(2,"Tsoil_K",   num_HRU,nYdirec_2openwq, nSnow_2openwq + nSoil_2openwq);

        // Master Json
        std::string master_json = std::getenv("master_json") ? std::getenv("master_json") : "";
        if (!std::filesystem::exists(master_json)) {
            std::cerr << "\nERROR: Path to OpenWQ_master.json does not exist !!\n"
                      << "Please set the environment variable 'master_json' "
                      << "to the path of the OpenWQ_master.json file.\n";
            exit(EXIT_FAILURE);
        }
        OpenWQ_wqconfig_ref->set_OpenWQ_masterjson(master_json);

        OpenWQ_couplercalls_ref->InitialConfig(
            *OpenWQ_hostModelconfig_ref,
            *OpenWQ_json_ref,                // create OpenWQ_json object
            *OpenWQ_wqconfig_ref,            // create OpenWQ_wqconfig object
            *OpenWQ_units_ref,               // functions for unit conversion
            *OpenWQ_utils_ref,                // utility methods/functions
            *OpenWQ_readjson_ref,            // read json files
            *OpenWQ_vars_ref,
            *OpenWQ_initiate_ref,            // initiate modules
            *OpenWQ_watertransp_ref,         // transport modules
            *OpenWQ_chem_ref,                // biochemistry modules
            *OpenWQ_extwatflux_ss_ref,       // sink and source modules)
            *OpenWQ_output_ref);
            
    }
    return 0;
}

// soilMoist_depVar does not have a value - it is passed as 0
int CLASSWQ_openwq::openwq_run_time_start(
    bool last_hru_flag,
    int index_hru, 
    int nSnow_2openwq, 
    int nSoil_2openwq,
    int simtime_summa[],
    double soilMoist_depVar_summa_frac[],                  
    double soilTemp_depVar_summa_K[],
    double airTemp_depVar_summa_K,
    double sweWatVol_stateVar_summa_m3[],
    double canopyWatVol_stateVar_summa_m3,
    double soilWatVol_stateVar_summa_m3[],
    double aquiferWatVol_stateVar_summa_m3) {

    time_t simtime = OpenWQ_units_ref->convertTime_ints2time_t(
        *OpenWQ_wqconfig_ref,
        simtime_summa[0], 
        simtime_summa[1], 
        simtime_summa[2], 
        simtime_summa[3], 
        simtime_summa[4],
        0);
    
    int runoff_vol = 0;
    
    // Updating Chemistry dependencies and volumes (out of order because of looping)

    OpenWQ_hostModelconfig_ref->set_dependVar_at(1,index_hru,0,0, airTemp_depVar_summa_K);
    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(canopy_index_openwq,index_hru,0,0, canopyWatVol_stateVar_summa_m3);   // canopy
    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(runoff_index_openwq,index_hru,0,0, runoff_vol);                       // runoff
    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(aquifer_index_openwq,index_hru,0,0, aquiferWatVol_stateVar_summa_m3); // aquifer

    // update Vars that rely on Snow
    for (int z = 0; z < nSnow_2openwq; z++) {
        OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(snow_index_openwq,index_hru,0,z, sweWatVol_stateVar_summa_m3[z]);  // snow
    }
    
    // Update Vars that rely on Soil
    for (int z = 0; z < nSoil_2openwq; z++) {
        OpenWQ_hostModelconfig_ref->set_dependVar_at(0,index_hru,0,z,soilMoist_depVar_summa_frac[z]); 
        OpenWQ_hostModelconfig_ref->set_dependVar_at(2,index_hru,0,z,soilTemp_depVar_summa_K[z]);
        OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(soil_index_openwq,index_hru,0,z, soilWatVol_stateVar_summa_m3[z]);      // soil

    }

    if (get_numHRU() -1 == index_hru ) {
        OpenWQ_couplercalls_ref->RunTimeLoopStart(
            *OpenWQ_hostModelconfig_ref,
            *OpenWQ_json_ref,
            *OpenWQ_wqconfig_ref,            // create OpenWQ_wqconfig object
            *OpenWQ_units_ref,               // functions for unit conversion
            *OpenWQ_utils_ref,                // utility methods/functions
            *OpenWQ_readjson_ref,            // read json files
            *OpenWQ_vars_ref,
            *OpenWQ_initiate_ref,            // initiate modules
            *OpenWQ_watertransp_ref,         // transport modules
            *OpenWQ_chem_ref,                // biochemistry modules
            *OpenWQ_extwatflux_ss_ref,          // sink and source modules)
            *OpenWQ_solver_ref,
            *OpenWQ_output_ref,
            simtime);
     }

    return 0;
}

int CLASSWQ_openwq::openwq_run_space(
    int simtime_summa[], 
    int source, int ix_s, int iy_s, int iz_s,
    int recipient, int ix_r, int iy_r, int iz_r, 
    double wflux_s2r, double wmass_source) {

    // Convert Fortran Index to C++ index
    ix_s -= 1; iy_s -= 1; iz_s -= 1;
    ix_r -= 1; iy_r -= 1; iz_r -= 1;

   
    time_t simtime = OpenWQ_units_ref->convertTime_ints2time_t(
        *OpenWQ_wqconfig_ref,
        simtime_summa[0], 
        simtime_summa[1], 
        simtime_summa[2], 
        simtime_summa[3], 
        simtime_summa[4],
        0);

    OpenWQ_couplercalls_ref->RunSpaceStep(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,            // create OpenWQ_wqconfig object
        *OpenWQ_units_ref,               // functions for unit conversion
        *OpenWQ_utils_ref,                // utility methods/functions
        *OpenWQ_readjson_ref,            // read json files
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,            // initiate modules
        *OpenWQ_watertransp_ref,         // transport modules
        *OpenWQ_chem_ref,                // biochemistry modules
        *OpenWQ_extwatflux_ss_ref,       // sink and source modules
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        simtime,
        source, ix_s, iy_s, iz_s,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r, wmass_source);

    return 0;
}

int CLASSWQ_openwq::openwq_run_space_in(
    int simtime_summa[],
    std::string source_EWF_name,
    int recipient, int ix_r, int iy_r, int iz_r, 
    double wflux_s2r) {

    // Convert Fortran Index to C++ index
    ix_r -= 1; iy_r -= 1; iz_r -= 1;
    
    time_t simtime = OpenWQ_units_ref->convertTime_ints2time_t(
        *OpenWQ_wqconfig_ref,
        simtime_summa[0], 
        simtime_summa[1], 
        simtime_summa[2], 
        simtime_summa[3], 
        simtime_summa[4],
        0);

     OpenWQ_couplercalls_ref->RunSpaceStep_IN(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,
        *OpenWQ_units_ref,
        *OpenWQ_utils_ref,
        *OpenWQ_readjson_ref,
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,
        *OpenWQ_watertransp_ref,
        *OpenWQ_chem_ref,
        *OpenWQ_extwatflux_ss_ref,
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        simtime,
        source_EWF_name,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r);

    return 0;
}

int CLASSWQ_openwq::openwq_run_time_end(
    int simtime_summa[]) {
    
    time_t simtime = OpenWQ_units_ref->convertTime_ints2time_t(
        *OpenWQ_wqconfig_ref,
        simtime_summa[0], 
        simtime_summa[1], 
        simtime_summa[2], 
        simtime_summa[3], 
        simtime_summa[4],
        0);


    OpenWQ_couplercalls_ref->RunTimeLoopEnd(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,            // create OpenWQ_wqconfig object
        *OpenWQ_units_ref,               // functions for unit conversion
        *OpenWQ_utils_ref,                // utility methods/functions
        *OpenWQ_readjson_ref,            // read json files
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,            // initiate modules
        *OpenWQ_watertransp_ref,         // transport modules
        *OpenWQ_chem_ref,                // biochemistry modules
        *OpenWQ_extwatflux_ss_ref,          // sink and source modules)
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        simtime);

    return 0;
}

int CLASSWQ_openwq::get_numHRU(){
    return this->num_HRU;
}
