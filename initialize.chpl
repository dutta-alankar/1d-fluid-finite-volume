module initialize {
    use globals;
    use structures;
    use Math;

    proc init_field (grid: borrowed Grid(?)): void {
        var xmin: (real(64), real(64), real(64)) = grid.xmin;
        var xmax: (real(64), real(64), real(64)) = grid.xmax;
        var extent: (real(64), real(64), real(64)) = (xmax[0]-xmin[0], xmax[1]-xmin[1], xmax[2]-xmin[2]);
        forall (i,j,k) in grid.indicesAll {
            var xc: (real(64), real(64), real(64)) = grid.cells_tot[i,j,k].center;
            var xl: (real(64), real(64), real(64)) = grid.cells_tot[i,j,k].wall_left.position;
            var xr: (real(64), real(64), real(64)) = grid.cells_tot[i,j,k].wall_right.position;
            // loop all the fluid variables
            for state_var in grid.cells_tot[i,j,k].state_prims_center.domain {
                if ((xc[0]-xmin[0])/extent[0])<(1.0/3.0) then
                    grid.cells_tot[i,j,k].state_prims_center[state_var] = 0.0;
                else if ((xc[0]-xmin[0])/extent[0])>=(1.0/3.0) && ((xc[0]-xmin[0])/extent[0])<(2.0/3.0) then
                    grid.cells_tot[i,j,k].state_prims_center[state_var] = 1.0;
                else
                    grid.cells_tot[i,j,k].state_prims_center[state_var] = 0.0;
            }
        }
        sync grid.cells_tot.updateFluff();
    }
}