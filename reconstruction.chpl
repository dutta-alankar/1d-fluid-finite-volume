module reconstruction {
    use globals;
    use structures;
    use utilities;
    use limiters;
    use Math;

    var reconstruction_coeff_domain: domain(1);
    var reconstruction_coeff: [reconstruction_coeff_domain] real(64);
    var use_slope_limiter: bool = true;
    
    proc compute_nghost(const ref reconstruction_type: string ): (int(64), int(64), int(64)) {
        if reconstruction_type=="constant" {
            reconstruction_coeff_domain = {1..1};
            for i in reconstruction_coeff_domain {
                reconstruction_coeff[i] = 1.0;
            }
        }
        else if reconstruction_type=="linear" {
            reconstruction_coeff_domain = {1..2};
            for i in reconstruction_coeff_domain {
                reconstruction_coeff[i] = 1.0; // i:real(64);
            }
        }
        else {
            writeln("Problem: Invaid reconstruction scheme!");
            exit(1);
        }
        return (reconstruction_coeff_domain.high, reconstruction_coeff_domain.high, reconstruction_coeff_domain.high);
    }

    proc interpolate_edges (grid: borrowed Grid(?)): void {
        sync grid.walls_tot.updateFluff();
        sync grid.cells_tot.updateFluff();
        var computeDomain: domain(3) = {grid.indicesInner.low[0]..grid.indicesInner.high[0], grid.indicesInner.low[1]..grid.indicesInner.high[1], grid.indicesInner.low[2]..grid.indicesInner.high[2]};
        // TODO: might need indicesInner
        forall (i,j,k) in grid.indicesAll {
            if i<(computeDomain.low[0]-1) || i>(computeDomain.high[0]+1) || j<(computeDomain.low[1]-1) || j>(computeDomain.high[1]+1) || k<(computeDomain.low[2]-1) || k>(computeDomain.high[2]+1) {
                for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                    grid.cells_tot[i,j,k].wall_left.state_consv_right[state_var] = 0.0;
                    grid.cells_tot[i,j,k].wall_right.state_consv_left[state_var] = 0.0;
                    grid.cells_tot[i,j,k].wall_front.state_consv_back[state_var] = 0.0;
                    grid.cells_tot[i,j,k].wall_back.state_consv_front[state_var] = 0.0;
                    grid.cells_tot[i,j,k].wall_bottom.state_consv_top[state_var] = 0.0;
                    grid.cells_tot[i,j,k].wall_top.state_consv_bottom[state_var] = 0.0;
                }
                continue;
            }
            // if !(grid.cells_tot[i,j,k].solve_flag && computeDomain.contains((i,j,k))) then continue;
            var sgn_dist_left: (real(64), real(64), real(64))  =(-(grid.cells_tot[i,j,k].center[0]-grid.cells_tot[i,j,k].wall_left.position[0])/grid.cells_tot[i,j,k].cell_size[0],
                                                                -(grid.cells_tot[i,j,k].center[1]-grid.cells_tot[i,j,k].wall_front.position[1])/grid.cells_tot[i,j,k].cell_size[1],
                                                                -(grid.cells_tot[i,j,k].center[2]-grid.cells_tot[i,j,k].wall_bottom.position[2])/grid.cells_tot[i,j,k].cell_size[2]); // negative
            var sgn_dist_right: (real(64), real(64), real(64)) = (-(grid.cells_tot[i,j,k].center[0]-grid.cells_tot[i,j,k].wall_right.position[0])/grid.cells_tot[i,j,k].cell_size[0],
                                                                -(grid.cells_tot[i,j,k].center[1]-grid.cells_tot[i,j,k].wall_back.position[1])/grid.cells_tot[i,j,k].cell_size[1],
                                                                -(grid.cells_tot[i,j,k].center[2]-grid.cells_tot[i,j,k].wall_top.position[2])/grid.cells_tot[i,j,k].cell_size[2]); // positive
            if sgn_dist_left[0]>0 || sgn_dist_left[1]>0 || sgn_dist_left[2]>0 then writeln("Problem on left distance");
            if sgn_dist_right[0]<0 || sgn_dist_right[1]<0 || sgn_dist_right[2]<0 then writeln("Problem on right distance");
            for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                var pow: int(64) = 0;
                var state_interp_left:  real(64) = 0.0;
                var state_interp_right: real(64) = 0.0;
                for l in reconstruction_coeff.domain {
                    var factor_x: real(64) = 0.0;
                    var factor_y: real(64) = 0.0;
                    var factor_z: real(64) = 0.0;
                    if l==reconstruction_coeff.domain.low then {
                        factor_x = grid.cells_tot[i,j,k].state_consv_center[state_var];
                        factor_y = grid.cells_tot[i,j,k].state_consv_center[state_var];
                        factor_z = grid.cells_tot[i,j,k].state_consv_center[state_var];
                    } else if l==(reconstruction_coeff.domain.low+1) {
                        if !use_slope_limiter then {
                            factor_x = 0.5*(grid.cells_tot[i+1,j,k].state_consv_center[state_var] - grid.cells_tot[i-1,j,k].state_consv_center[state_var]);
                            factor_y = 0.5*(grid.cells_tot[i,j+1,k].state_consv_center[state_var] - grid.cells_tot[i,j-1,k].state_consv_center[state_var]);
                            factor_z = 0.5*(grid.cells_tot[i,j,k+1].state_consv_center[state_var] - grid.cells_tot[i,j,k-1].state_consv_center[state_var]);
                        } else if slope_limiter=="minmod" then {
                            factor_x = minmod(grid.cells_tot[i,j,k].state_consv_center[state_var] - grid.cells_tot[i-1,j,k].state_consv_center[state_var],
                                            grid.cells_tot[i+1,j,k].state_consv_center[state_var] - grid.cells_tot[i,j,k].state_consv_center[state_var]);
                            factor_y = minmod(grid.cells_tot[i,j,k].state_consv_center[state_var] - grid.cells_tot[i,j-1,k].state_consv_center[state_var],
                                            grid.cells_tot[i,j+1,k].state_consv_center[state_var] - grid.cells_tot[i,j,k].state_consv_center[state_var]);
                            factor_z = minmod(grid.cells_tot[i,j,k].state_consv_center[state_var] - grid.cells_tot[i,j,k-1].state_consv_center[state_var],
                                            grid.cells_tot[i,j,k+1].state_consv_center[state_var] - grid.cells_tot[i,j,k].state_consv_center[state_var]);
                        } else if slope_limiter=="mc" then {
                            factor_x = mc_lim(0.5*(grid.cells_tot[i+1,j,k].state_consv_center[state_var] - grid.cells_tot[i-1,j,k].state_consv_center[state_var]),
                                            2*(grid.cells_tot[i+1,j,k].state_consv_center[state_var] - grid.cells_tot[i,j,k].state_consv_center[state_var]),
                                            2*(grid.cells_tot[i,j,k].state_consv_center[state_var] - grid.cells_tot[i-1,j,k].state_consv_center[state_var]));
                            factor_y = mc_lim(0.5*(grid.cells_tot[i,j+1,k].state_consv_center[state_var] - grid.cells_tot[i,j-1,k].state_consv_center[state_var]),
                                            2*(grid.cells_tot[i,j+1,k].state_consv_center[state_var] - grid.cells_tot[i,j,k].state_consv_center[state_var]),
                                            2*(grid.cells_tot[i,j,k].state_consv_center[state_var] - grid.cells_tot[i,j-1,k].state_consv_center[state_var]));
                            factor_z = mc_lim(0.5*(grid.cells_tot[i,j,k+1].state_consv_center[state_var] - grid.cells_tot[i,j,k-1].state_consv_center[state_var]),
                                            2*(grid.cells_tot[i,j,k+1].state_consv_center[state_var] - grid.cells_tot[i,j,k].state_consv_center[state_var]),
                                            2*(grid.cells_tot[i,j,k].state_consv_center[state_var] - grid.cells_tot[i,j,k-1].state_consv_center[state_var]));
                        } else {
                            writeln("Slope limiter ", slope_limiter, " is not supported!");
                            exit(1);
                        }
                    } else {
                        factor_x = 0.0;
                        factor_y = 0.0;
                        factor_z = 0.0;
                    }
                    // if debug then writeln(grid.cells_tot[i,j,k].center, " l=", l, ", pow=", pow, ": ", factor_x * reconstruction_coeff[l] * sgn_dist_left[0]**pow); 
                    state_interp_left  +=  (factor_x * reconstruction_coeff[l] * sgn_dist_left[0]**pow);
                    state_interp_right +=  (factor_x * reconstruction_coeff[l] * sgn_dist_right[0]**pow);
                    pow += 1;
                }
                grid.cells_tot[i,j,k].wall_left.state_consv_right[state_var] = state_interp_left;
                grid.cells_tot[i,j,k].wall_right.state_consv_left[state_var] = state_interp_right;
                grid.cells_tot[i,j,k].wall_front.state_consv_back[state_var] = state_interp_left; // Placeholder, needs proper calculation
                grid.cells_tot[i,j,k].wall_back.state_consv_front[state_var] = state_interp_right; // Placeholder, needs proper calculation
                grid.cells_tot[i,j,k].wall_bottom.state_consv_top[state_var] = state_interp_left; // Placeholder, needs proper calculation
                grid.cells_tot[i,j,k].wall_top.state_consv_bottom[state_var] = state_interp_right; // Placeholder, needs proper calculation
            }
        }
        // sync grid.cells_tot.updateFluff();
        // sync grid.walls_tot.updateFluff();
        
        forall (i,j,k) in grid.indicesAll {
            if i<computeDomain.low[0] || i>(computeDomain.high[0]+1) || j<computeDomain.low[1] || j>(computeDomain.high[1]+1) || k<computeDomain.low[2] || k>(computeDomain.high[2]+1) {
                for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                    grid.walls_tot[i,j,k].state_consv_left[state_var]  = 0.0;
                    grid.walls_tot[i,j,k].state_consv_right[state_var]  = 0.0;
                }
                continue;
            }
            for state_var in grid.cells_tot[i,j,k].state_consv_center.domain {
                grid.walls_tot[i,j,k].state_consv_left[state_var]  = grid.cells_tot[i-1,j,k].wall_right.state_consv_left[state_var];
                grid.walls_tot[i,j,k].state_consv_right[state_var] = grid.cells_tot[i,j,k].wall_left.state_consv_right[state_var];
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