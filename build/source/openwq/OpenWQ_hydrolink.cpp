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

// OpenWQ hydrolink for SUMMA with (optional) internally coupled mizuRoute.
// See OpenWQ_hydrolink.h for the compartments.

#include "OpenWQ_hydrolink.h"
#include "OpenWQ_interface.h"

#include <algorithm>
#include <fstream>

CLASSWQ_openwq::CLASSWQ_openwq() {}

CLASSWQ_openwq::~CLASSWQ_openwq() {}

// SUMMA time (year, month, day, hour, minute) to time_t
time_t CLASSWQ_openwq::to_time(int simtime_summa[]) {
    return OpenWQ_units_ref->convertTime_ints2time_t(
        *OpenWQ_wqconfig_ref,
        simtime_summa[0],
        simtime_summa[1],
        simtime_summa[2],
        simtime_summa[3],
        simtime_summa[4],
        0);
}

// =============================================================================
// decl: declare compartments, external fluxes, exports and dependencies
// =============================================================================
// colId[x]  : id of the HRU of land column x
// colDom[x] : 0 if the HRU has a single domain, otherwise the domain index (1-based)
// reachId[] : ids of the river reaches (num_reach = 0 without a river network)
int CLASSWQ_openwq::decl(
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

    this->num_col   = num_col;
    this->num_reach = num_reach;

    if (OpenWQ_hostModelconfig_ref->get_num_HydroComp() != 0) return 0;

    // -------------------------------------------------------------------------
    // Compartments (capital letters)
    // -------------------------------------------------------------------------
    std::vector<int> nx_cmp, nz_cmp;

    OpenWQ_hostModelconfig_ref->add_HydroComp(
        canopy_index_openwq, "SCALARCANOPYWAT",
        num_col, nYdirec_2openwq, nCanopy_2openwq);
    nx_cmp.push_back(num_col); nz_cmp.push_back(nCanopy_2openwq);

    OpenWQ_hostModelconfig_ref->add_HydroComp(
        snow_index_openwq, "ILAYERVOLFRACWAT_SNOW",
        num_col, nYdirec_2openwq, nSnow_2openwq);
    nx_cmp.push_back(num_col); nz_cmp.push_back(nSnow_2openwq);

    OpenWQ_hostModelconfig_ref->add_HydroComp(
        runoff_index_openwq, "RUNOFF",
        num_col, nYdirec_2openwq, nRunoff_2openwq);
    nx_cmp.push_back(num_col); nz_cmp.push_back(nRunoff_2openwq);

    OpenWQ_hostModelconfig_ref->add_HydroComp(
        soil_index_openwq, "ILAYERVOLFRACWAT_SOIL",
        num_col, nYdirec_2openwq, nSoil_2openwq);
    nx_cmp.push_back(num_col); nz_cmp.push_back(nSoil_2openwq);

    OpenWQ_hostModelconfig_ref->add_HydroComp(
        aquifer_index_openwq, "SCALARAQUIFER",
        num_col, nYdirec_2openwq, nAquifer_2openwq);
    nx_cmp.push_back(num_col); nz_cmp.push_back(nAquifer_2openwq);

    // water on its way to the stream: runoff delivered by the columns, held while SUMMA routes it
    OpenWQ_hostModelconfig_ref->add_HydroComp(
        stream_index_openwq, "RUNOFF_TO_STREAM",
        num_col, nYdirec_2openwq, 1);
    nx_cmp.push_back(num_col); nz_cmp.push_back(1);

    int next_index = stream_index_openwq + 1;

    if (nLake_2openwq > 0) {
        lake_index_openwq = next_index++;
        OpenWQ_hostModelconfig_ref->add_HydroComp(
            lake_index_openwq, "ILAYERVOLFRACWAT_LAKE",
            num_col, nYdirec_2openwq, nLake_2openwq);
        nx_cmp.push_back(num_col); nz_cmp.push_back(nLake_2openwq);
    }

    if (num_reach > 0) {
        river_index_openwq = next_index++;
        OpenWQ_hostModelconfig_ref->add_HydroComp(
            river_index_openwq, "RIVER_NETWORK_REACHES",
            num_reach, 1, 1);
        nx_cmp.push_back(num_reach); nz_cmp.push_back(1);
    }

    // -------------------------------------------------------------------------
    // External water fluxes (capital letters)
    // -------------------------------------------------------------------------
    OpenWQ_hostModelconfig_ref->add_HydroExtFlux(
        0, "PRECIP", num_col, nYdirec_2openwq, 1);
    if (has_glacier != 0) {
        OpenWQ_hostModelconfig_ref->add_HydroExtFlux(
            1, "GLACIER_ICE_MELT", num_col, nYdirec_2openwq, 1);
    }

    // -------------------------------------------------------------------------
    // Flux-concentration exports (selected with FLUXES_CONC_TO_PRINT in the master file).
    // Each land export is named after the SUMMA variable that gives its volume.
    // -------------------------------------------------------------------------
    OpenWQ_hostModelconfig_ref->add_FluxConcExport(
        scalarRunoffVol_fluxexp_openwq, "scalarRunoffVol_m3",
        runoff_index_openwq, num_col, nYdirec_2openwq, nRunoff_2openwq);
    OpenWQ_hostModelconfig_ref->add_FluxConcExport(
        averageRoutedRunoff_fluxexp_openwq, "averageRoutedRunoff",
        stream_index_openwq, num_col, nYdirec_2openwq, 1);
    OpenWQ_hostModelconfig_ref->add_FluxConcExport(
        scalarTotalRunoff_fluxexp_openwq, "scalarTotalRunoff",
        stream_index_openwq, num_col, nYdirec_2openwq, 1);
    if (num_reach > 0) {
        OpenWQ_hostModelconfig_ref->add_FluxConcExport(
            reachOutflow_fluxexp_openwq, "Qlocal_out",
            river_index_openwq, num_reach, 1, 1);
    }

    OpenWQ_vars_ref = std::make_unique<OpenWQ_vars>(
        OpenWQ_hostModelconfig_ref->get_num_HydroComp(),
        OpenWQ_hostModelconfig_ref->get_num_HydroExtFlux());

    // -------------------------------------------------------------------------
    // Dependency variables. They are indexed by cell (ix, iy, iz) whatever the
    // compartment, so all share one size that covers every compartment.
    // A reach shares its values with the land column of the same index when
    // there is one (see openwq_set_reach_state).
    // -------------------------------------------------------------------------
    const int nx_depend = std::max(num_col, num_reach);
    nz_depend = std::max({nSnow_2openwq + nSoil_2openwq, nLake_2openwq, 1});

    OpenWQ_hostModelconfig_ref->add_HydroDepend(
        sm_depend_openwq, "SM", nx_depend, nYdirec_2openwq, nz_depend);
    OpenWQ_hostModelconfig_ref->add_HydroDepend(
        tair_depend_openwq, "Tair_K", nx_depend, nYdirec_2openwq, nz_depend);
    OpenWQ_hostModelconfig_ref->add_HydroDepend(
        tsoil_depend_openwq, "Tsoil_K", nx_depend, nYdirec_2openwq, nz_depend);
    OpenWQ_hostModelconfig_ref->add_HydroDepend(
        swrad_depend_openwq, "SWrad_Wm2", nx_depend, nYdirec_2openwq, nz_depend);
    OpenWQ_hostModelconfig_ref->add_HydroDepend(
        area_depend_openwq, "cellArea_m2", nx_depend, nYdirec_2openwq, nz_depend);
    OpenWQ_hostModelconfig_ref->add_HydroDepend(
        tc_depend_openwq, "T", nx_depend, nYdirec_2openwq, nz_depend);

    if (num_reach > 0 && num_col > 1) {
        std::cout << "<OpenWQ> WARNING: the dependency variables (temperature, radiation, area) of the first "
                  << std::min(num_col, num_reach)
                  << " river reaches are those of the land column with the same index." << std::endl;
    }

    // -------------------------------------------------------------------------
    // Cell ids, so that the SS/EWF files and the output refer to SUMMA and
    // mizuRoute ids. Must be set before InitialConfig().
    //   land:  "<hruId>_z<layer>", or "<hruId>_d<domain>_z<layer>" when the HRU has several domains
    //   river: "<segId>"
    // -------------------------------------------------------------------------
    OpenWQ_hostModelconfig_ref->set_cellid_to_wqlabel("hruId");

    for (int cmp = 0; cmp < (int)OpenWQ_hostModelconfig_ref->get_num_HydroComp(); cmp++) {

        arma::Cube<double> domain_xyz(nx_cmp[cmp], nYdirec_2openwq, nz_cmp[cmp]);
        OpenWQ_hostModelconfig_ref->set_cellid_to_wq_size(domain_xyz);

        for (int x = 0; x < nx_cmp[cmp]; x++) {
            if (cmp == river_index_openwq) {
                OpenWQ_hostModelconfig_ref->set_cellid_to_wq_at(
                    cmp, x, 0, 0,
                    std::to_string(static_cast<long long>(reachId[x])));
                continue;
            }
            std::string base = std::to_string(static_cast<long long>(colId[x]));
            if (colDom[x] > 0) base += "_d" + std::to_string(colDom[x]);
            for (int z = 0; z < nz_cmp[cmp]; z++) {
                OpenWQ_hostModelconfig_ref->set_cellid_to_wq_at(
                    cmp, x, 0, z,
                    base + "_z" + std::to_string(z + 1));
            }
        }
    }

    OpenWQ_wqconfig_ref->set_OpenWQ_masterjson("openWQ_master.json");

    OpenWQ_couplercalls_ref->InitialConfig(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,
        *OpenWQ_units_ref,
        *OpenWQ_utils_ref,
        *OpenWQ_readjson_ref,
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,
        *OpenWQ_TD_model_ref,
        *OpenWQ_LE_model_ref,
        *OpenWQ_CH_model_ref,
        *OpenWQ_SI_model_ref,
        *OpenWQ_TS_model_ref,
        *OpenWQ_extwatflux_ss_ref,
        *OpenWQ_output_ref);

    // Compartment index table for the supporting scripts: spatial parameters
    // address compartments by index, and land and river cells share index ranges
    try {
        std::string out_dir = OpenWQ_wqconfig_ref->get_output_dir();
        if (!out_dir.empty()) {
            std::filesystem::create_directories(out_dir);
            std::ofstream f(out_dir + "/openwq_compartments.json");
            int ncmp = OpenWQ_hostModelconfig_ref->get_num_HydroComp();
            f << "{\n  \"host\": \"summa" << (num_reach > 0 ? "+mizuroute" : "") << "\",\n  \"compartments\": [\n";
            for (int c = 0; c < ncmp; c++) {
                f << "    {\"index\": " << c << ", \"name\": \""
                  << OpenWQ_hostModelconfig_ref->get_HydroComp_name_at(c)
                  << "\", \"nx\": " << nx_cmp[c] << ", \"nz\": " << nz_cmp[c] << "}"
                  << (c + 1 < ncmp ? "," : "") << "\n";
            }
            f << "  ],\n  \"river_compartment\": " << river_index_openwq
              << ",\n  \"lake_compartment\": " << lake_index_openwq << "\n}\n";
        }
    } catch (...) {}

    // after InitialConfig (memory allocated) and after the cell ids are registered
    OpenWQ_couplercalls_ref->ParseEWFandSS(
        *OpenWQ_json_ref,
        *OpenWQ_vars_ref,
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_wqconfig_ref,
        *OpenWQ_units_ref,
        *OpenWQ_utils_ref,
        *OpenWQ_output_ref,
        *OpenWQ_extwatflux_ss_ref);

    return 0;
}

// =============================================================================
// openwq_set_reach_state: reach water volume and dependencies at the start of the step
// =============================================================================
// Call before the land columns: a reach whose index is also a land column keeps
// the dependency values of that column.
int CLASSWQ_openwq::openwq_set_reach_state(
    int n_reach,
    double reachVol_m3[],
    double airTemp_K[],
    double SWrad_Wm2[],
    double area_m2[]) {

    if (river_index_openwq < 0) return 0;

    for (int r = 0; r < n_reach; r++) {

        OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
            river_index_openwq, r, 0, 0, reachVol_m3[r]);

        if (r < num_col) continue;

        for (int z = 0; z < nz_depend; z++) {
            OpenWQ_hostModelconfig_ref->set_dependVar_at(sm_depend_openwq,    r, 0, z, 1.0);
            OpenWQ_hostModelconfig_ref->set_dependVar_at(tair_depend_openwq,  r, 0, z, airTemp_K[r]);
            OpenWQ_hostModelconfig_ref->set_dependVar_at(tsoil_depend_openwq, r, 0, z, airTemp_K[r]);
            OpenWQ_hostModelconfig_ref->set_dependVar_at(swrad_depend_openwq, r, 0, z, SWrad_Wm2[r]);
            OpenWQ_hostModelconfig_ref->set_dependVar_at(area_depend_openwq,  r, 0, z, area_m2[r]);
            OpenWQ_hostModelconfig_ref->set_dependVar_at(tc_depend_openwq,    r, 0, z, airTemp_K[r] - 273.15);
        }
    }

    return 0;
}

// =============================================================================
// openwq_run_time_start: water volumes and dependencies of one land column
// =============================================================================
// index_col is 0-based. The last column starts the OpenWQ time step.
int CLASSWQ_openwq::openwq_run_time_start(
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

    // dependencies: every layer gets a value (below the soil, the deepest soil layer)
    for (int z = 0; z < nz_depend; z++) {
        const int zs = std::min(z, nSoil_2openwq - 1);
        const double tsoil = (zs >= 0) ? soilTemp_depVar_summa_K[zs] : airTemp_depVar_summa_K;
        const double sm    = (zs >= 0) ? soilMoist_depVar_summa_frac[zs] : 0.0;
        OpenWQ_hostModelconfig_ref->set_dependVar_at(sm_depend_openwq,    index_col, 0, z, sm);
        OpenWQ_hostModelconfig_ref->set_dependVar_at(tair_depend_openwq,  index_col, 0, z, airTemp_depVar_summa_K);
        OpenWQ_hostModelconfig_ref->set_dependVar_at(tsoil_depend_openwq, index_col, 0, z, tsoil);
        OpenWQ_hostModelconfig_ref->set_dependVar_at(swrad_depend_openwq, index_col, 0, z, SWrad_depVar_summa_Wm2);
        OpenWQ_hostModelconfig_ref->set_dependVar_at(area_depend_openwq,  index_col, 0, z, col_area_m2);
        OpenWQ_hostModelconfig_ref->set_dependVar_at(tc_depend_openwq,    index_col, 0, z, tsoil - 273.15);
    }

    // single-layer compartments; the two pools hold no water until the fluxes are known
    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
        canopy_index_openwq, index_col, 0, 0, canopyWatVol_stateVar_summa_m3);
    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
        runoff_index_openwq, index_col, 0, 0, 0.0);
    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
        aquifer_index_openwq, index_col, 0, 0, aquiferWatVol_stateVar_summa_m3);
    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
        stream_index_openwq, index_col, 0, 0, 0.0);

    // layered compartments; layers that do not exist now hold no water
    {
        const int nz = (int) OpenWQ_hostModelconfig_ref->get_HydroComp_num_cells_z_at(snow_index_openwq);
        for (int z = 0; z < nz; z++) {
            const double vol = (z < nSnow_2openwq) ? sweWatVol_stateVar_summa_m3[z] : 0.0;
            OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
                snow_index_openwq, index_col, 0, z, vol);
        }
    }
    {
        const int nz = (int) OpenWQ_hostModelconfig_ref->get_HydroComp_num_cells_z_at(soil_index_openwq);
        for (int z = 0; z < nz; z++) {
            const double vol = (z < nSoil_2openwq) ? soilWatVol_stateVar_summa_m3[z] : 0.0;
            OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
                soil_index_openwq, index_col, 0, z, vol);
        }
    }
    if (lake_index_openwq >= 0) {
        const int nz = (int) OpenWQ_hostModelconfig_ref->get_HydroComp_num_cells_z_at(lake_index_openwq);
        for (int z = 0; z < nz; z++) {
            const double vol = (z < nLake_2openwq) ? lakeWatVol_stateVar_summa_m3[z] : 0.0;
            OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
                lake_index_openwq, index_col, 0, z, vol);
        }
    }

    // start the time step once all columns are updated
    if (index_col == num_col - 1) {
        OpenWQ_couplercalls_ref->RunTimeLoopStart(
            *OpenWQ_hostModelconfig_ref,
            *OpenWQ_json_ref,
            *OpenWQ_wqconfig_ref,
            *OpenWQ_units_ref,
            *OpenWQ_utils_ref,
            *OpenWQ_readjson_ref,
            *OpenWQ_vars_ref,
            *OpenWQ_initiate_ref,
            *OpenWQ_TD_model_ref,
            *OpenWQ_LE_model_ref,
            *OpenWQ_CH_model_ref,
            *OpenWQ_SI_model_ref,
            *OpenWQ_TS_model_ref,
            *OpenWQ_extwatflux_ss_ref,
            *OpenWQ_solver_ref,
            *OpenWQ_output_ref,
            to_time(simtime_summa));
    }

    return 0;
}

// =============================================================================
// openwq_run_space: water flux between two cells
// =============================================================================
int CLASSWQ_openwq::openwq_run_space(
    int simtime_summa[],
    int source, int ix_s, int iy_s, int iz_s,
    int recipient, int ix_r, int iy_r, int iz_r,
    double wflux_s2r, double wmass_source,
    int out_to_stream) {

    // 1-based to 0-based; -1 (out of the domain) is kept
    ix_s = std::max(-1, ix_s - 1);
    iy_s = std::max(-1, iy_s - 1);
    iz_s = std::max(-1, iz_s - 1);
    ix_r = std::max(-1, ix_r - 1);
    iy_r = std::max(-1, iy_r - 1);
    iz_r = std::max(-1, iz_r - 1);

    // Water that SUMMA sends to the stream leaves the land column (recipient -1) with the
    // transport of the source compartment unchanged; the dissolved mass it carries is then
    // placed in the RUNOFF_TO_STREAM cell of the same column instead of leaving the domain.
    const bool to_stream = (out_to_stream != 0)
        && (recipient == -1)
        && (source != stream_index_openwq)
        && (source != river_index_openwq)
        && (ix_s >= 0) && (iy_s >= 0) && (iz_s >= 0);

    std::vector<double> d_source_before;
    if (to_stream) {
        const auto& d_src = (*OpenWQ_vars_ref->d_chemass_dt_transp_diss)(source);
        d_source_before.resize(d_src.n_elem);
        for (unsigned int chemi = 0; chemi < d_src.n_elem; chemi++)
            d_source_before[chemi] = d_src(chemi)(ix_s, iy_s, iz_s);
    }

    OpenWQ_couplercalls_ref->RunSpaceStep(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,
        *OpenWQ_units_ref,
        *OpenWQ_utils_ref,
        *OpenWQ_readjson_ref,
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,
        *OpenWQ_TD_model_ref,
        *OpenWQ_TS_model_ref,
        *OpenWQ_LE_model_ref,
        *OpenWQ_CH_model_ref,
        *OpenWQ_SI_model_ref,
        *OpenWQ_extwatflux_ss_ref,
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        to_time(simtime_summa),
        source, ix_s, iy_s, iz_s,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r, wmass_source);

    if (to_stream) {
        const auto& d_src = (*OpenWQ_vars_ref->d_chemass_dt_transp_diss)(source);
        auto& d_pool = (*OpenWQ_vars_ref->d_chemass_dt_transp_diss)(stream_index_openwq);
        auto& mb = OpenWQ_vars_ref->mass_balance;
        for (unsigned int chemi = 0; chemi < d_src.n_elem; chemi++) {
            const double delivered = d_source_before[chemi] - d_src(chemi)(ix_s, iy_s, iz_s);
            if (delivered <= 0.0) continue;
            d_pool(chemi)(ix_s, 0, 0) += delivered;
            // still in the domain: it is counted as an outflow when it leaves the pool or the river
            if (mb.initialized && chemi < mb.num_species)
                mb.cumulative_out_flux[chemi] -= delivered;
        }
    }

    return 0;
}

// =============================================================================
// openwq_run_space_in: water flux entering the domain (e.g. precipitation)
// =============================================================================
int CLASSWQ_openwq::openwq_run_space_in(
    int simtime_summa[],
    std::string source_EWF_name,
    int recipient, int ix_r, int iy_r, int iz_r,
    double wflux_s2r) {

    ix_r -= 1;
    iy_r -= 1;
    iz_r -= 1;

    OpenWQ_couplercalls_ref->RunSpaceStep_IN(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,
        *OpenWQ_units_ref,
        *OpenWQ_utils_ref,
        *OpenWQ_readjson_ref,
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,
        *OpenWQ_TD_model_ref,
        *OpenWQ_CH_model_ref,
        *OpenWQ_TS_model_ref,
        *OpenWQ_extwatflux_ss_ref,
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        to_time(simtime_summa),
        source_EWF_name,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r);

    return 0;
}

// =============================================================================
// openwq_set_watervol: water volume of one cell
// =============================================================================
// Used for the through-flow pools and the reaches, whose mixing volume is only
// known once the fluxes of the step are known.
int CLASSWQ_openwq::openwq_set_watervol(
    int icmp, int ix, int iy, int iz, double vol_m3) {

    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
        icmp, ix - 1, iy - 1, iz - 1, vol_m3);

    return 0;
}

// =============================================================================
// openwq_set_fluxvol: through-volume of a flux-concentration export
// =============================================================================
// The output writer forms conc = mass/volume of the source cell and mass = conc * this volume.
int CLASSWQ_openwq::openwq_set_fluxvol(
    int iflux, int ix, int iy, int iz, double flux_vol_m3) {

    OpenWQ_hostModelconfig_ref->set_fluxVol_hydromodel_at(
        iflux, ix - 1, iy - 1, iz - 1, flux_vol_m3);

    return 0;
}

// =============================================================================
// openwq_run_time_end: solve the time step and write the output
// =============================================================================
int CLASSWQ_openwq::openwq_run_time_end(
    int simtime_summa[]) {

    OpenWQ_couplercalls_ref->RunTimeLoopEnd(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,
        *OpenWQ_units_ref,
        *OpenWQ_utils_ref,
        *OpenWQ_readjson_ref,
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,
        *OpenWQ_TD_model_ref,
        *OpenWQ_LE_model_ref,
        *OpenWQ_CH_model_ref,
        *OpenWQ_SI_model_ref,
        *OpenWQ_TS_model_ref,
        *OpenWQ_extwatflux_ss_ref,
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        to_time(simtime_summa));

    return 0;
}
