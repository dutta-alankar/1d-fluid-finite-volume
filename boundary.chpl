module boundary {
    use globals;
    use structures;
    use conversion;

    proc set_boundary (grid: borrowed Grid(?)): void {
        var computeDomain: domain(3) = {grid.indicesInner.low[0]..grid.indicesInner.high[0], grid.indicesInner.low[1]..grid.indicesInner.high[1], grid.indicesInner.low[2]..grid.indicesInner.high[2]};
        forall (i,j,k) in grid.indicesAll {
            if computeDomain.contains((i,j,k)) {
                /* place code here to set internal boundaries */
                /*
                (i,j,k) th cell needs to be set some primitive value as an user-defined boundary
                for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                        grid.cells_tot[i,j,k].state_prims_center[state_var] =  <set what you want>; 
                }
                prims_to_consv(grid, only_cell=true, cell_id=(i,j,k))
                */
            }
            else {
                // X-direction boundaries
                if i < computeDomain.low[0] {
                    for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                        grid.cells_tot[i,j,k].state_prims_center[state_var] =  grid.cells_tot[i+grid.npoints_int[0],j,k].state_prims_center[state_var]; // periodic
                    }
                } else if i > computeDomain.high[0] {
                    for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                        grid.cells_tot[i,j,k].state_prims_center[state_var] =  grid.cells_tot[i-grid.npoints_int[0],j,k].state_prims_center[state_var]; // periodic
                    }
                }
                // Y-direction boundaries
                if j < computeDomain.low[1] {
                    for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                        grid.cells_tot[i,j,k].state_prims_center[state_var] =  grid.cells_tot[i,j+grid.npoints_int[1],k].state_prims_center[state_var]; // periodic
                    }
                } else if j > computeDomain.high[1] {
                    for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                        grid.cells_tot[i,j,k].state_prims_center[state_var] =  grid.cells_tot[i,j-grid.npoints_int[1],k].state_prims_center[state_var]; // periodic
                    }
                }
                // Z-direction boundaries
                if k < computeDomain.low[2] {
                    for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                        grid.cells_tot[i,j,k].state_prims_center[state_var] =  grid.cells_tot[i,j,k+grid.npoints_int[2]].state_prims_center[state_var]; // periodic
                    }
                } else if k > computeDomain.high[2] {
                    for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                        grid.cells_tot[i,j,k].state_prims_center[state_var] =  grid.cells_tot[i,j,k-grid.npoints_int[2]].state_prims_center[state_var]; // periodic
                    }
                }
                prims_to_consv(grid, only_boundary=true);
            }
        }

        sync grid.cells_tot.updateFluff();
        sync grid.walls_tot.updateFluff();
    }
}