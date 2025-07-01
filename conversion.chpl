module conversion {
    use globals;
    use structures;

    proc consv_to_prims (grid: borrowed Grid(?), only_boundary: bool=false, only_cell: bool=false, cell_id:(int(64), int(64), int(64))=(0,0,0)): void {
        var computeDomain: domain(3) = {grid.indicesInner.low[0]..grid.indicesInner.high[0], grid.indicesInner.low[1]..grid.indicesInner.high[1], grid.indicesInner.low[2]..grid.indicesInner.high[2]};
        // TODO: might need indicesInner
        forall (i,j,k) in grid.indicesAll {
            if only_boundary && computeDomain.contains((i,j,k)) then continue;
            if only_cell && (i,j,k)!=cell_id then continue;
            for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                grid.cells_tot[i,j,k].state_prims_center[state_var] =  grid.cells_tot[i,j,k].state_consv_center[state_var];
            }
        }
        if !only_cell then sync grid.cells_tot.updateFluff();
    }

    proc prims_to_consv (grid: borrowed Grid(?), only_boundary: bool=false, only_cell: bool=false, cell_id:(int(64), int(64), int(64))=(0,0,0)): void {
        var computeDomain: domain(3) = {grid.indicesInner.low[0]..grid.indicesInner.high[0], grid.indicesInner.low[1]..grid.indicesInner.high[1], grid.indicesInner.low[2]..grid.indicesInner.high[2]};
        // TODO: might need indicesInner
        forall (i,j,k) in grid.indicesAll {
            if only_boundary && computeDomain.contains((i,j,k)) then continue;
            if only_cell && (i,j,k)!=cell_id then continue;
            for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                grid.cells_tot[i,j,k].state_consv_center[state_var] =  grid.cells_tot[i,j,k].state_prims_center[state_var];
            }
        }
        if !only_cell then sync grid.cells_tot.updateFluff();
    }
}