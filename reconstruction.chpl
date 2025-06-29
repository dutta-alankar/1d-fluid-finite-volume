module reconstruction {
    use globals;
    use structures;
    use Math;

    var reconstruction_coeff_domain: domain(1);
    var reconstruction_coeff: [reconstruction_coeff_domain] real(64);
    
    proc compute_nghost(const ref reconstruction_type: string ): int(64) {
        if reconstruction_type=="constant" {
            reconstruction_coeff_domain = {1..1};
            for i in reconstruction_coeff_domain {
                reconstruction_coeff[i] = 1.0;
            }
        }
        else if reconstruction_type=="linear" {
            reconstruction_coeff_domain = {1..2};
            for i in reconstruction_coeff_domain {
                reconstruction_coeff[i] = 1.0/i:real(64);
            }
        }
        else {
            writeln("Problem: Invaid reconstruction scheme!");
            exit(1);
        }
        return reconstruction_coeff_domain.high;
    }

    proc interpolate_edges (grid: borrowed Grid(?)): void {
        sync grid.walls_tot.updateFluff();
        sync grid.cells_tot.updateFluff();
        var computeDomain: domain(1) = {grid.indicesInner.low..grid.indicesInner.high};
        // TODO: might need indicesInner
        forall i in grid.indicesAll {
            if i<(computeDomain.low-1) || i>(computeDomain.high+1) {
                for state_var in grid.cells_tot[i].state_consv_center.domain {
                    grid.cells_tot[i].wall_left.state_consv_right[state_var] = 0.0;
                    grid.cells_tot[i].wall_right.state_consv_left[state_var] = 0.0;
                }
                continue;
            }
            // if !(grid.cells_tot[i].solve_flag && computeDomain.contains(i)) then continue;
            var sgn_dist_left: real(64)  = -(grid.cells_tot[i].center-grid.cells_tot[i].wall_left.position)/grid.cells_tot[i].cell_size; // negative
            var sgn_dist_right: real(64) = -(grid.cells_tot[i].center-grid.cells_tot[i].wall_right.position)/grid.cells_tot[i].cell_size; // positive
            if sgn_dist_left>0 then writeln("Problem on left distance");
            if sgn_dist_right<0 then writeln("Problem on right distance");
            for state_var in grid.cells_tot[i].state_consv_center.domain {
                var pow: int(64) = 0;
                var state_interp_left:  real(64) = 0.0;
                var state_interp_right: real(64) = 0.0;
                for j in reconstruction_coeff.domain {
                    var factor: real(64) = 0.0;
                    if j==reconstruction_coeff.domain.low then 
                        factor = grid.cells_tot[i].state_consv_center[state_var];
                    else if j==(reconstruction_coeff.domain.low+1) then {
                        factor = grid.cells_tot[i+1].state_consv_center[state_var] - grid.cells_tot[i-1].state_consv_center[state_var];
                    }
                    else
                        factor = 0.0;
                    if debug then writeln(grid.cells_tot[i].center, " j=", j, ", pow=", pow, ": ", factor * reconstruction_coeff[j] * sgn_dist_left**pow); 
                    state_interp_left  +=  (factor * reconstruction_coeff[j] * sgn_dist_left**pow);
                    state_interp_right +=  (factor * reconstruction_coeff[j] * sgn_dist_right**pow);
                    pow += 1;
                }
                grid.cells_tot[i].wall_left.state_consv_right[state_var] = state_interp_left;
                grid.cells_tot[i].wall_right.state_consv_left[state_var] = state_interp_right;
            }
        }
        // sync grid.cells_tot.updateFluff();
        // sync grid.walls_tot.updateFluff();
        
        forall i in grid.indicesAll {
            if i<computeDomain.low || i>(computeDomain.high+1) {
                for state_var in grid.cells_tot[i].state_consv_center.domain {
                    grid.walls_tot[i].state_consv_left[state_var]  = 0.0;
                    grid.walls_tot[i].state_consv_right[state_var]  = 0.0;
                }
                continue;
            }
            for state_var in grid.cells_tot[i].state_consv_center.domain {
                grid.walls_tot[i].state_consv_left[state_var]  = grid.cells_tot[i-1].wall_right.state_consv_left[state_var];
                grid.walls_tot[i].state_consv_right[state_var] = grid.cells_tot[i].wall_left.state_consv_right[state_var];
            }
        }
        sync grid.walls_tot.updateFluff();
        sync grid.cells_tot.updateFluff();
        /*
        for i in grid.indicesAll {
            if i<computeDomain.low || i>(computeDomain.high+1) then continue;
            for state_var in grid.cells_tot[i].state_consv_center.domain {
                write("wall: i=", i, ", x=", grid.walls_tot[i].position, ": ", 
                        grid.walls_tot[i].state_consv_left[state_var], ", ", 
                        grid.walls_totindicesAll[i].state_consv_right[state_var]);
                if i<=computeDomain.high then
                    writeln(" cell ", grid.cells_tot[i-1].center, ": ", 
                            grid.cells_tot[i-1].state_consv_center[state_var], ", ",
                            grid.cells_tot[i].center, ": ",
                            grid.cells_tot[i].state_consv_center[state_var]);
                else
                    writeln(" cell ", grid.cells_tot[i-1].center, ": ", grid.cells_tot[i-1].state_consv_center[state_var]);
            }
        }*/
    }
}