module initialize {
    use structures;

    proc init_field (grid: borrowed Grid(?)): void {
        var xmin: real(64) = grid.xmin;
        var xmax: real(64) = grid.xmax;
        var extent: real(64) = xmax-xmin;
        sync forall i in grid.indicesInner {
            var xc: real(64) = grid.cells_tot[i].center;
            var xl: real(64) = grid.cells_tot[i].wall_left.position;
            var xr: real(64) = grid.cells_tot[i].wall_right.position;
            // loop all the fluid variables
            for nvar in grid.cells_tot[i].state_prims_center.domain {
                if ((xc-xmin)/extent)<=0.5 then
                    grid.cells_tot[i].state_prims_center[nvar] = 1.0;
                else
                    grid.cells_tot[i].state_prims_center[nvar] = 0.0;
            }
        }
        sync grid.cells_tot.updateFluff();
    }
}