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

// C interface of the OpenWQ hydrolink (called from openWQ.f90 through iso_c_binding).
// The arguments are described in OpenWQ_hydrolink.h.

#ifndef OPENWQ_INTERFACE_H
#define OPENWQ_INTERFACE_H

#ifdef __cplusplus
extern "C" {
    class CLASSWQ_openwq;
    typedef CLASSWQ_openwq CLASSWQ_openwq;
#else
    typedef struct CLASSWQ_openwq CLASSWQ_openwq;
#endif

    CLASSWQ_openwq* create_openwq();

    void delete_openwq(CLASSWQ_openwq* openWQ);

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
        long long reachId[]);

    int openwq_set_reach_state(
        CLASSWQ_openwq *openWQ,
        int n_reach,
        double reachVol_m3[],
        double airTemp_K[],
        double SWrad_Wm2[],
        double area_m2[]);

    int openwq_run_time_start(
        CLASSWQ_openwq *openWQ,
        int index_col,
        int nSnow_2openwq,
        int nLake_2openwq,
        int nSoil_2openwq,
        int simtime_summa[],
        double soilMoist_depVar[],
        double soilTemp_K_depVar[],
        double airTemp_K_depVar,
        double SWrad_Wm2_depVar,
        double sweWatVol_stateVar[],
        double lakeWatVol_stateVar[],
        double canopyWat,
        double soilWatVol_stateVar[],
        double aquiferStorage,
        double col_area_m2);

    int openwq_run_space(
        CLASSWQ_openwq *openWQ,
        int simtime_summa[],
        int source, int ix_s, int iy_s, int iz_s,
        int recipient, int ix_r, int iy_r, int iz_r,
        double wflux_s2r, double wmass_source,
        int out_to_stream);

    int openwq_run_space_in(
        CLASSWQ_openwq *openWQ,
        int simtime_summa[],
        char* source_EWF_name,
        int recipient, int ix_r, int iy_r, int iz_r,
        double wflux_s2r);

    int openwq_set_watervol(
        CLASSWQ_openwq *openWQ,
        int icmp, int ix, int iy, int iz,
        double vol_m3);

    int openwq_set_fluxvol(
        CLASSWQ_openwq *openWQ,
        int iflux, int ix, int iy, int iz,
        double flux_vol_m3);

    int openwq_run_time_end(
        CLASSWQ_openwq *openWQ,
        int simtime_summa[]);

#ifdef __cplusplus
}
#endif

#endif // OPENWQ_INTERFACE_H
