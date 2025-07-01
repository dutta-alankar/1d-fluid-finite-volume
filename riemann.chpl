module riemann {
    use globals;
    use structures;

    /* riemann solver for the advection problem with constant velocity */
    proc solver_advect (const ref wall: Wall, vel: real(64), direction: int(64)): [state_domain] real(64) {
        var solution_consv: [state_domain] real(64);
        for state_var in state_domain {
            if vel>=0 then {
                select direction {
                    when 0 do solution_consv[state_var] = wall.state_consv_left[state_var];
                    when 1 do solution_consv[state_var] = wall.state_consv_front[state_var];
                    when 2 do solution_consv[state_var] = wall.state_consv_bottom[state_var];
                    otherwise do halt("Invalid direction");
                }
            } else {
               select direction {
                    when 0 do solution_consv[state_var] = wall.state_consv_right[state_var];
                    when 1 do solution_consv[state_var] = wall.state_consv_back[state_var];
                    when 2 do solution_consv[state_var] = wall.state_consv_top[state_var];
                    otherwise do halt("Invalid direction");
                }
            }
        }
        return solution_consv; 
    }

    proc calc_flux (solution_consv: [state_domain] real(64), vel: real(64)): [state_domain] real(64) {
        var flux_solve: [state_domain] real(64);
        for state_var in state_domain {
            flux_solve[state_var] = vel*solution_consv[state_var];
        }
        return flux_solve;
    }

    proc solve_at_walls (grid: borrowed Grid(?), vel: (real(64), real(64), real(64))): void {
        var computeDomain: domain(3) = {grid.indicesInner.low[0]..grid.indicesInner.high[0], grid.indicesInner.low[1]..grid.indicesInner.high[1], grid.indicesInner.low[2]..grid.indicesInner.high[2]};
        forall (i,j,k) in grid.indicesAll {
            if !grid.cells_tot[i,j,k].solve_flag && !(i>=computeDomain.low[0] && i<=(computeDomain.high[0]+1) && j>=computeDomain.low[1] && j<=(computeDomain.high[1]+1) && k>=computeDomain.low[2] && k<=(computeDomain.high[2]+1)) then {
                continue;
            }

            // Solve in X-direction
            var solution_consv_x: [state_domain] real(64) = solver_advect (grid.walls_tot[i,j,k], vel[0], 0);
            var flux_solve_x: [state_domain] real(64) = calc_flux(solution_consv_x, vel[0]);
            for state_var in state_domain {
                grid.walls_tot[i,j,k].state_consv_solve[state_var] =  solution_consv_x[state_var];
                grid.walls_tot[i,j,k].flux_solve[state_var]  = flux_solve_x[state_var];
                if grid.cells_tot[i,j,k].solve_flag {
                    grid.cells_tot[i,j,k].wall_left.state_consv_solve[state_var] =  solution_consv_x[state_var];
                    grid.cells_tot[i,j,k].wall_left.flux_solve[state_var]  = flux_solve_x[state_var];
                }
                if grid.cells_tot[i-1,j,k].solve_flag {
                    grid.cells_tot[i-1,j,k].wall_right.state_consv_solve[state_var] =  solution_consv_x[state_var];
                    grid.cells_tot[i-1,j,k].wall_right.flux_solve[state_var]  = flux_solve_x[state_var];
                }
            }

            // Solve in Y-direction
            var solution_consv_y: [state_domain] real(64) = solver_advect (grid.walls_tot[i,j,k], vel[1], 1);
            var flux_solve_y: [state_domain] real(64) = calc_flux(solution_consv_y, vel[1]);
            for state_var in state_domain {
                if grid.cells_tot[i,j,k].solve_flag {
                    grid.cells_tot[i,j,k].wall_front.state_consv_solve[state_var] =  solution_consv_y[state_var];
                    grid.cells_tot[i,j,k].wall_front.flux_solve[state_var]  = flux_solve_y[state_var];
                }
                if grid.cells_tot[i,j-1,k].solve_flag {
                    grid.cells_tot[i,j-1,k].wall_back.state_consv_solve[state_var] =  solution_consv_y[state_var];
                    grid.cells_tot[i,j-1,k].wall_back.flux_solve[state_var]  = flux_solve_y[state_var];
                }
            }

            // Solve in Z-direction
            var solution_consv_z: [state_domain] real(64) = solver_advect (grid.walls_tot[i,j,k], vel[2], 2);
            var flux_solve_z: [state_domain] real(64) = calc_flux(solution_consv_z, vel[2]);
            for state_var in state_domain {
                if grid.cells_tot[i,j,k].solve_flag {
                    grid.cells_tot[i,j,k].wall_bottom.state_consv_solve[state_var] =  solution_consv_z[state_var];
                    grid.cells_tot[i,j,k].wall_bottom.flux_solve[state_var]  = flux_solve_z[state_var];
                }
                if grid.cells_tot[i,j,k-1].solve_flag {
                    grid.cells_tot[i,j,k-1].wall_top.state_consv_solve[state_var] =  solution_consv_z[state_var];
                    grid.cells_tot[i,j,k-1].wall_top.flux_solve[state_var]  = flux_solve_z[state_var];
                }
            }
        }
        sync grid.walls_tot.updateFluff();
        sync grid.cells_tot.updateFluff();
    }
}