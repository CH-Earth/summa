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

// C wrappers of the CLASSWQ_openwq methods, callable from Fortran.

#include "OpenWQ_hydrolink.h"
#include "OpenWQ_interface.h"

CLASSWQ_openwq* create_openwq() {
    return new CLASSWQ_openwq();
}

void delete_openwq(CLASSWQ_openwq* openWQ) {
    delete openWQ;
}

int openwq_decl(
    CLASSWQ_openwq *openWQ,
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
    long long reachId[]) {

    return openWQ->decl(
        num_col,
        nCanopy_2openwq,
        nSnow_2openwq,
        nSoil_2openwq,
        nRunoff_2openwq,
        nAquifer_2openwq,
        nLake_2openwq,
        nYdirec_2openwq,
        colId,
        colDom,
        has_glacier,
        num_reach,
        reachId);
}

int openwq_set_reach_state(
    CLASSWQ_openwq *openWQ,
    int n_reach,
    double reachVol_m3[],
    double airTemp_K[],
    double SWrad_Wm2[],
    double area_m2[]) {

    return openWQ->openwq_set_reach_state(
        n_reach, reachVol_m3, airTemp_K, SWrad_Wm2, area_m2);
}

int openwq_run_time_start(
    CLASSWQ_openwq *openWQ,
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
    double col_area_m2) {

    return openWQ->openwq_run_time_start(
        index_col,
        nSnow_2openwq,
        nLake_2openwq,
        nSoil_2openwq,
        simtime_summa,
        soilMoist_depVar_summa_frac,
        soilTemp_depVar_summa_K,
        airTemp_depVar_summa_K,
        SWrad_depVar_summa_Wm2,
        sweWatVol_stateVar_summa_m3,
        lakeWatVol_stateVar_summa_m3,
        canopyWatVol_stateVar_summa_m3,
        soilWatVol_stateVar_summa_m3,
        aquiferWatVol_stateVar_summa_m3,
        col_area_m2);
}

int openwq_run_space(
    CLASSWQ_openwq *openWQ,
    int simtime_summa[],
    int source, int ix_s, int iy_s, int iz_s,
    int recipient, int ix_r, int iy_r, int iz_r,
    double wflux_s2r, double wmass_source,
    int out_to_stream) {

    return openWQ->openwq_run_space(
        simtime_summa,
        source, ix_s, iy_s, iz_s,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r, wmass_source,
        out_to_stream);
}

int openwq_run_space_in(
    CLASSWQ_openwq *openWQ,
    int simtime_summa[],
    char* source_EWF_name,
    int recipient, int ix_r, int iy_r, int iz_r,
    double wflux_s2r) {

    std::string source_EWF_name_str(source_EWF_name);

    return openWQ->openwq_run_space_in(
        simtime_summa,
        source_EWF_name_str,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r);
}

int openwq_set_watervol(
    CLASSWQ_openwq *openWQ,
    int icmp, int ix, int iy, int iz,
    double vol_m3) {

    return openWQ->openwq_set_watervol(icmp, ix, iy, iz, vol_m3);
}

int openwq_set_fluxvol(
    CLASSWQ_openwq *openWQ,
    int iflux, int ix, int iy, int iz,
    double flux_vol_m3) {

    return openWQ->openwq_set_fluxvol(iflux, ix, iy, iz, flux_vol_m3);
}

int openwq_run_time_end(
    CLASSWQ_openwq *openWQ,
    int simtime_summa[]) {

    return openWQ->openwq_run_time_end(simtime_summa);
}
