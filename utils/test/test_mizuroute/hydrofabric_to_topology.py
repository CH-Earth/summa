#!/usr/bin/env python3
"""Build a mizuRoute topology file from a real NextGen hydrofabric.

make_test_topology.py exercises the SUMMA-mizuRoute coupling with a synthetic,
single-chain network that carries no hydrological meaning. This script instead
reads the real river network out of a NextGen hydrofabric geopackage (the same
one t-route reads directly, e.g. gage-10154200_subset.gpkg for the Provo test
domain) and writes it in the minimal topology format the SUMMA-mizuRoute
coupling understands (see make_test_topology.py / config.toml [hydrofabric]).

That makes it possible to route the same SUMMA runoff through the same real
network with mizuRoute and with t-route (via ngen), and compare the two --
see compare_to_troute.py and the folder README.

Hydrofabric layout (NextGen "hydrofabric" v2 gpkg):
  - flowpaths / flowpath-attributes: one reach per catchment, id 'wb-<n>',
    draining to a nexus 'nex-<m>'.
  - nexus: 'nex-<m>' draining to the next flowpath downstream (or to a
    flowpath outside the subset, which marks the domain outlet).
  - divides: one catchment 'cat-<n>' per flowpath, with drainage area.

SUMMA/NextGen catchment <-> GRU mapping:
  Each catchment's cat-<n>.input file sets attrib_file_HRU_order, the
  1-based position of that catchment's GRU in the domain's attributes.nc
  (SUMMA and the hydrofabric use unrelated id numbering, so this index is
  the only link between them).

Usage:
    hydrofabric_to_topology.py <domain_dir> <topology.nc>

<domain_dir> is a test_ngen domain folder (e.g. domain_provo) containing
settings/<gage>.gpkg, settings/SUMMA/attributes.nc and settings/cat-*.input.
"""

import glob
import os
import re
import sqlite3
import sys

import numpy as np
from netCDF4 import Dataset

HRU_ORDER_RE = re.compile(r"attrib_file_HRU_order\s*=\s*(\d+)")
CAT_ID_RE = re.compile(r"cat-(\d+)\.input$")


def find_gpkg(settings_dir):
    matches = glob.glob(os.path.join(settings_dir, "*.gpkg"))
    if len(matches) != 1:
        sys.exit(f"expected exactly one .gpkg in {settings_dir}, found {matches}")
    return matches[0]


def read_hru_order(settings_dir):
    """catchment number -> 1-based index into attributes.nc arrays."""
    order = {}
    for f in glob.glob(os.path.join(settings_dir, "cat-*.input")):
        cat_id = int(CAT_ID_RE.search(f).group(1))
        m = HRU_ORDER_RE.search(open(f).read())
        if not m:
            sys.exit(f"no attrib_file_HRU_order in {f}")
        order[cat_id] = int(m.group(1))
    return order


def read_network(gpkg_file):
    """Return {catchment_id: (down_catchment_id_or_0, length_m, slope, area_m2)}."""
    con = sqlite3.connect(gpkg_file)
    cur = con.cursor()

    cur.execute("SELECT id, toid FROM flowpaths")
    fp_toid = dict(cur.fetchall())  # 'wb-<n>' -> 'nex-<m>'

    cur.execute("SELECT id, toid FROM nexus")
    nex_toid = dict(cur.fetchall())  # 'nex-<m>' -> 'wb-<n>' (or out of subset)

    cur.execute('SELECT id, "Length_m", "So" FROM "flowpath-attributes"')
    fp_attrib = {wb: (length, slope) for wb, length, slope in cur.fetchall()}

    cur.execute("SELECT divide_id, areasqkm FROM divides")
    divide_area = {int(cid.split("-")[1]): area for cid, area in cur.fetchall()}

    con.close()

    network = {}
    n_outlets = 0
    for wb, nex in fp_toid.items():
        cat_id = int(wb.split("-")[1])

        down_wb = nex_toid.get(nex)
        if down_wb is not None and down_wb in fp_toid:
            down_id = int(down_wb.split("-")[1])
        else:
            down_id = 0
            n_outlets += 1

        length, slope = fp_attrib[wb]
        area = divide_area[cat_id] * 1.0e6  # km2 -> m2

        network[cat_id] = (down_id, length, max(slope, 1.0e-5), area)

    if n_outlets != 1:
        sys.exit(f"expected exactly one outlet reach, found {n_outlets}")

    return network


def main(domain_dir, topology_file):
    settings_dir = os.path.join(domain_dir, "settings")
    summa_dir = os.path.join(settings_dir, "SUMMA")
    gpkg_file = find_gpkg(settings_dir)

    print(f"hydrofabric:  {gpkg_file}")

    hru_order = read_hru_order(settings_dir)
    network = read_network(gpkg_file)

    if set(hru_order) != set(network):
        only_summa = set(hru_order) - set(network)
        only_hf = set(network) - set(hru_order)
        sys.exit(
            "catchment id mismatch between cat-*.input and the hydrofabric\n"
            f"  in cat-*.input only: {sorted(only_summa)}\n"
            f"  in hydrofabric only: {sorted(only_hf)}"
        )

    n = len(hru_order)

    with Dataset(os.path.join(summa_dir, "attributes.nc")) as attrib:
        gru_id = np.array(attrib.variables["gruId"][:], dtype="i8")
        hru_area_summa = np.array(attrib.variables["HRUarea"][:], dtype="f8")

    if len(gru_id) != n:
        sys.exit(f"attributes.nc has {len(gru_id)} GRUs, hydrofabric has {n} catchments")

    cat_ids = sorted(hru_order)
    hru_id = np.array([gru_id[hru_order[c] - 1] for c in cat_ids], dtype="i8")
    seg_id = np.array(cat_ids, dtype="i8")
    hru_to_seg_id = seg_id.copy()  # each catchment's runoff enters its own reach
    down_seg_id = np.array([network[c][0] for c in cat_ids], dtype="i8")
    length = np.array([network[c][1] for c in cat_ids], dtype="f8")
    slope = np.array([network[c][2] for c in cat_ids], dtype="f8")
    area = np.array([network[c][3] for c in cat_ids], dtype="f8")

    # sanity check: hydrofabric divide area should match the SUMMA HRU area
    # (independent confirmation that attrib_file_HRU_order lines the two up correctly)
    summa_area = np.array([hru_area_summa[hru_order[c] - 1] for c in cat_ids])
    rel_err = np.abs(area - summa_area) / summa_area
    if np.any(rel_err > 0.05):
        bad = [(c, a, s) for c, a, s, r in zip(cat_ids, area, summa_area, rel_err) if r > 0.05]
        print("WARNING: hydrofabric and SUMMA areas disagree by >5% for:")
        for c, a, s in bad:
            print(f"  cat-{c}: hydrofabric={a:.4e} m2, SUMMA={s:.4e} m2")
    else:
        print(f"area cross-check: hydrofabric and SUMMA agree to within {rel_err.max()*100:.2f}%")

    with Dataset(topology_file, "w", format="NETCDF4") as dst:
        dst.createDimension("hru", n)
        dst.createDimension("seg", n)

        def put(name, dim, data, dtype, units, long_name):
            var = dst.createVariable(name, dtype, (dim,))
            var[:] = data
            var.units = units
            var.long_name = long_name

        put("hruId", "hru", hru_id, "i8", "-", "routing HRU id, matched to SUMMA GRU ids")
        put("area", "hru", area, "f8", "m2", "routing HRU (catchment) area, from the hydrofabric")
        put("hruToSegId", "hru", hru_to_seg_id, "i8", "-", "id of the stream segment below each HRU")
        put("segId", "seg", seg_id, "i8", "-", "stream segment id (hydrofabric catchment/flowpath number)")
        put("downSegId", "seg", down_seg_id, "i8", "-", "downstream segment id (0 at the outlet)")
        put("length", "seg", length, "f8", "m", "segment length, from flowpath-attributes.Length_m")
        put("slope", "seg", slope, "f8", "-", "segment slope, from flowpath-attributes.So")

        dst.description = (
            "real river network for the NextGen hydrofabric domain, converted for the "
            "SUMMA-mizuRoute coupling by hydrofabric_to_topology.py. segId/hruToSegId "
            "equal the hydrofabric catchment number (e.g. segId=2863621 is wb-2863621 / "
            "cat-2863621), so mizuRoute's Q_reach can be lined up with t-route output "
            "for the same hydrofabric."
        )
        dst.source_hydrofabric = os.path.abspath(gpkg_file)

    print(f"wrote {topology_file}: {n} routing HRUs, {n} reaches")
    outlet = cat_ids[int(np.where(down_seg_id == 0)[0][0])]
    print(f"outlet: cat-{outlet} / wb-{outlet}")
    print(f"total area: {area.sum():.4e} m2")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
