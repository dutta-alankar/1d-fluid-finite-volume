module riemann {
    use structures;

    /* riemann solver for the advection problem with constant velocity */
    proc solver_advect (const ref wall: Wall, vel: real(64)): [state_domain] real(64) {
        var solution_consv: [state_domain] real(64);
        for state_var in state_domain {
            if vel>=0 then
                solution_consv[state_var] = wall.state_consv_left[state_var];
            else
               solution_consv[state_var] = wall.state_consv_right[state_var]; 
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

    proc solve_at_walls (grid: borrowed Grid(?), vel: real(64)): void {
        var computeDomain: domain(1) = {grid.indicesInner.low..grid.indicesInner.high};
        sync forall i in grid.indicesAll {
            // if debug then writeln("i = ", i, ", flag = ", grid.cells_tot[i].solve_flag, " ", computeDomain.high+1);
            if !grid.cells_tot[i].solve_flag && i!=(computeDomain.high+1) then continue;
            var solution_consv: [state_domain] real(64) = solver_advect (grid.walls_tot[i], vel);
            var flux_solve: [state_domain] real(64) = calc_flux(solution_consv, vel);
            /*
            for state_var in state_domain {
                if solution_consv[state_var]<0 then solution_consv[state_var] = 0.0;
                if flux_solve[state_var]<0 then flux_solve[state_var] = 0.0;
            }
            */
            if debug then writeln("i = ", i, ": consv_state: ", solution_consv, " flux: ", flux_solve);
            for state_var in state_domain {
                grid.walls_tot[i].state_consv_solve[state_var] =  solution_consv[state_var];
                grid.walls_tot[i].flux_solve[state_var]  = flux_solve[state_var];
                if grid.cells_tot[i].solve_flag {
                    grid.cells_tot[i].wall_left.state_consv_solve[state_var] =  solution_consv[state_var];
                    grid.cells_tot[i].wall_left.flux_solve[state_var]  = flux_solve[state_var];
                }
                if grid.cells_tot[i-1].solve_flag {
                    grid.cells_tot[i-1].wall_right.state_consv_solve[state_var] =  solution_consv[state_var];
                    grid.cells_tot[i-1].wall_right.flux_solve[state_var]  = flux_solve[state_var];
                }
            }
        }
        sync grid.walls_tot.updateFluff();
        sync grid.cells_tot.updateFluff();
    }
}