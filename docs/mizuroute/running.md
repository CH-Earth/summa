# Running mizuRoute

mizuRoute can be run with SUMMA using either the **sequential** or **coupled**
configuration described in [Coupling design](design.md).

## Sequential configuration

In the sequential configuration, SUMMA and mizuRoute are run as separate
executables. SUMMA is run first to generate runoff, and the resulting SUMMA
NetCDF output is subsequently used as input to standalone mizuRoute.

Configuration of a standalone mizuRoute simulation is described in the
[mizuRoute documentation](https://mizuroute.readthedocs.io/en/main/users_guide/).
The standalone simulation is controlled using the standard mizuRoute control
file.

For example, a control file for a SUMMA--mizuRoute application specifies the
simulation period and routing method, the river-network topology, and the
SUMMA NetCDF file and variable containing the runoff supplied to mizuRoute.
It also specifies whether runoff must be remapped between the SUMMA spatial
units and the river-network catchments.

Once the mizuRoute control file has been prepared, the standalone model can be
run with:

```bash
mpirun -np <n> route_runoff <mizuRoute_control_file>
```

where `<n>` is the number of MPI processes, `route_runoff` is the standalone
mizuRoute executable, and `<mizuRoute_control_file>` is the mizuRoute control
file.

For example:

```bash
mpirun -np 4 route_runoff settings/mizuRoute/mizuroute.control
```

Users should refer to the
[mizuRoute documentation](https://mizuroute.readthedocs.io/en/main/users_guide/)
for details of the standalone control file, model options, and parallel
execution.

## Coupled configuration

In the coupled configuration, mizuRoute is executed directly within the SUMMA
simulation. SUMMA is run normally, with an additional TOML configuration file
provided using the `-c` or `--config` command-line option.

For example:

```bash
summa.exe \
    -m <summa_control_file> \
    -c <summa_config_file>
```

The existing SUMMA control file specified with `-m` is unchanged. The
additional configuration file specified with `-c` currently contains the
information required to initialize and run mizuRoute.

The TOML configuration is currently used primarily to provide the additional
information required by mizuRoute. The longer-term plan is to migrate SUMMA
control information to TOML, as described in
[Coupling architecture](coupling.md).

An example configuration is:

```toml
# mizuRoute control file

[mizuRoute]

namelist_path = "settings/mizuRoute/"
namelist_file = "param.nml.default"

methods = "3"

dt = 3600

[hydrofabric]

hfabric_path = "settings/mizuRoute/"
hfabric_file = "topology.nc"

dname_seg = "seg"
dname_hru = "hru"

varname_HRUid     = "hruId"
varname_segId     = "segId"
varname_hruSegId  = "hruToSegId"
varname_downSegId = "downSegId"

varname_area   = "area"
varname_slope  = "slope"
varname_length = "length"
```

The `[mizuRoute]` section defines the routing configuration. It specifies the
mizuRoute parameter namelist, the routing method or methods to execute, and the
routing time step.

The `[hydrofabric]` section identifies the river-network topology file and
defines the dimension and variable names used to describe the river network,
including:

- river reaches and hydrologic response units;
- HRU and reach identifiers;
- the reach associated with each HRU;
- downstream reach connectivity; and
- catchment area, channel slope, and channel length.

The `methods` option identifies the mizuRoute routing method to use. In the
example above, `methods = "3"` selects the Eulerian kinematic-wave routing
method. See the
[mizuRoute documentation](https://mizuroute.readthedocs.io/)
for descriptions of the available routing methods.

The `dt` option specifies the mizuRoute routing time step in seconds. In the
example above, `dt = 3600` specifies an hourly routing time step. The routing 
time step must be less than or equal to the SUMMA model time step.

### Coupled output

mizuRoute output is written directly into the SUMMA NetCDF output files. The
coupled configuration adds a `seg` dimension for river reaches and a `method`
dimension for the routing method. Routed streamflow is written as:

```text
q_reach(time, seg, method)
```

where q\_reach is streamflow at the downstream end of each river reach, with
units of m3 s-1. The seg coordinate contains the river-reach identifiers,
and upArea(seg) provides the drainage area above the downstream end of each
reach.

Unlike the sequential configuration, the coupled configuration does not
require a mizuRoute runoff-input file. SUMMA passes the simulated runoff
directly to mizuRoute during model execution. Routed streamflow is
written into the SUMMA NetCDF output rather than to a separate standalone
mizuRoute output file.

## Parallel execution

The current coupled SUMMA--mizuRoute configuration does not support parallel
execution of the coupled land--river system. Users who require mizuRoute's
parallel routing capabilities should use the **sequential configuration** and
run the standalone mizuRoute executable using MPI.

The reasons for this current limitation, and the intended approach for
parallelization of the coupled system, are described in
[Coupling design](design.md#parallel-execution).
