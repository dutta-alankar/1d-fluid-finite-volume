module boundary {
    use structures;
    use conversion;

    proc set_boundary (grid: borrowed Grid(?)): void {
        var computeDomain: domain(1) = {grid.indicesInner.low..grid.indicesInner.high};
        sync forall i in grid.indicesAll {
            if computeDomain.contains(i) {
                /* place code here to set internal boundaries */
                /*
                i th cell needs to be set some primitive value as an user-defined boundary
                for state_var in grid.cells_tot[i].state_consv_center.domain {
                        grid.cells_tot[i].state_prims_center[state_var] =  <set what you want>; 
                }
                prims_to_consv(grid, only_cell=true, cell_id=i)
                */
            }
            else {
                if i<computeDomain.low {
                    /* place code here to set the left boundary */
                    for state_var in grid.cells_tot[i].state_consv_center.domain {
                        grid.cells_tot[i].state_prims_center[state_var] =  grid.cells_tot[i+grid.npoints_int].state_prims_center[state_var]; // periodic
                    }
                    // if debug then writeln("L ", i, ": ", grid.cells_tot[i].state_prims_center, " <- R ", i+grid.npoints_int, ": ", grid.cells_tot[i+grid.npoints_int].state_prims_center);
                }
                else {
                    /* place code here to set the right boundary */
                    for state_var in grid.cells_tot[i].state_consv_center.domain {
                        grid.cells_tot[i].state_prims_center[state_var] =  grid.cells_tot[i%grid.npoints_int].state_prims_center[state_var]; // periodic
                    }
                    // if debug then writeln("R ", i, ": ", grid.cells_tot[i].state_prims_center, " <- L ", i%grid.npoints_int, ": ", grid.cells_tot[i%grid.npoints_int].state_prims_center);
                }
                prims_to_consv(grid, only_boundary=true);
            }
        }

        sync grid.cells_tot.updateFluff();
    }
}