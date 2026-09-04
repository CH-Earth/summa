# mizuRoute

[mizuRoute](https://github.com/ESCOMP/mizuRoute) is a river network routing model
designed to route runoff from hydrologic and land surface models through
large-scale river networks. It supports vector-based river networks and
provides multiple river-routing methods for simulating streamflow throughout
the river network.

SUMMA can be built with mizuRoute to route runoff directly during a SUMMA
simulation. This provides an integrated workflow in which SUMMA simulates the
terrestrial water balance and mizuRoute routes the resulting runoff through
the river network.

These pages describe how to build and run SUMMA with mizuRoute. For details
on mizuRoute itself, including its routing methods, input data, and model
configuration, see the
[mizuRoute documentation](https://mizuroute.readthedocs.io/).
