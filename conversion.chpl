module conversion {
    use structures;

    proc consv_to_prims (grid: borrowed Grid(?), only_boundary: bool=false, only_cell: bool=false, cell_id:int(64)=0): void {
        var computeDomain: domain(1) = {grid.indicesInner.low..grid.indicesInner.high};
        // TODO: might need indicesInner
        forall i in grid.indicesAll {
            if only_boundary && computeDomain.contains(i) then continue;
            if only_cell && i!=cell_id then continue;
            for state_var in grid.cells_tot[i].state_consv_center.domain {
                grid.cells_tot[i].state_prims_center[state_var] =  grid.cells_tot[i].state_consv_center[state_var];
            }
        }
        if !only_cell then sync grid.cells_tot.updateFluff();
    }

    proc prims_to_consv (grid: borrowed Grid(?), only_boundary: bool=false, only_cell: bool=false, cell_id:int(64)=0): void {
        var computeDomain: domain(1) = {grid.indicesInner.low..grid.indicesInner.high};
        // TODO: might need indicesInner
        forall i in grid.indicesAll {
            if only_boundary && computeDomain.contains(i) then continue;
            if only_cell && i!=cell_id then continue;
            for state_var in grid.cells_tot[i].state_consv_center.domain {
                grid.cells_tot[i].state_consv_center[state_var] =  grid.cells_tot[i].state_prims_center[state_var];
            }
        }
        if !only_cell then sync grid.cells_tot.updateFluff();
    }
}