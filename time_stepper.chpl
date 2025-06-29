module time_stepper {
    use globals;
    use structures;

    proc timeStepEuler (grid: borrowed Grid(?), dt: real(64)): void {
        /* explicit Euler first order time integration */
        // riemann problem solve
        solve_at_walls(grid, vel_adv);
        forall i in grid.indicesInner {
            if !grid.cells_tot[i].solve_flag then continue;
            for state_var in grid.cells_tot[i].states_count {
                grid.cells_tot[i].rhs_state[state_var] = (grid.cells_tot[i].wall_left.flux_solve[state_var] - grid.cells_tot[i].wall_right.flux_solve[state_var])/grid.cells_tot[i].cell_size;
            }
        }
        forall i in grid.indicesInner {
            if !grid.cells_tot[i].solve_flag then continue;
            for state_var in grid.cells_tot[i].states_count {
                grid.cells_tot[i].state_consv_center[state_var] = grid.cells_tot[i].state_consv_center[state_var] + grid.cells_tot[i].rhs_state[state_var]*dt;
            }
        }
    }

    proc timeStepRK2 (grid: borrowed Grid(?), dt: real(64)): void {
        /* explicit RK2 second order time integration */
        var grid_tmp: Grid(?) = new Grid();
        // RK2 step 1 : use the already calculated flux for half timeStep
        writeln("RK2 step 1");
        grid.deepCopy(grid_tmp);
        timeStepEuler(grid_tmp, 0.5*dt);
        consv_to_prims(grid_tmp);

        write("k1: [");
        for i in grid_tmp.indicesInner {
            if !grid_tmp.cells_tot[i].solve_flag then continue;
            for state_var in grid_tmp.cells_tot[i].states_count {
                write(grid_tmp.cells_tot[i].rhs_state[state_var], ", ");
            }
        }
        writeln("\b\b]");
        // set the rhs fluxes based on half timestep value of conservative states at half timestep
        /*
        forall i in grid_tmp.indicesInner {
            if !grid_tmp.cells_tot[i].solve_flag then continue;
            for state_var in grid_tmp.cells_tot[i].states_count {
                grid_tmp.cells_tot[i].state_consv_center[state_var] = grid_tmp.cells_tot[i].state_consv_center[state_var] + grid_tmp.cells_tot[i].rhs_state[state_var]*dt;
            }
        }
        sync grid_tmp.cells_tot.updateFluff();
        sync grid_tmp.walls_tot.updateFluff();
        */
        // RK2 step 2 : use the half timeStep conservative states to evaluate conservative states
        writeln("RK2 step 2:");
        set_boundary(grid_tmp);
        interpolate_edges(grid_tmp);
        // riemann problem solve
        solve_at_walls(grid_tmp, vel_adv);
        forall i in grid.indicesInner {
            if !grid.cells_tot[i].solve_flag then continue;
            for state_var in grid.cells_tot[i].states_count {
                grid.cells_tot[i].rhs_state[state_var] = (grid_tmp.cells_tot[i].wall_left.flux_solve[state_var] - grid_tmp.cells_tot[i].wall_right.flux_solve[state_var])/grid_tmp.cells_tot[i].cell_size;
            }
        }
        write("k2: [");
        for i in grid.indicesInner {
            if !grid.cells_tot[i].solve_flag then continue;
            for state_var in grid.cells_tot[i].states_count {
                write(grid.cells_tot[i].rhs_state[state_var], ", ");
            }
        }
        writeln("\b\b]");
        
        // riemann problem solve
        // solve_at_walls(grid_tmp, vel_adv);
        // set the rhs fluxes based on half timestep value of conservative states
        forall i in grid.indicesInner {
            if !grid.cells_tot[i].solve_flag then continue;
            for state_var in grid.cells_tot[i].states_count {
                grid.cells_tot[i].state_consv_center[state_var] = grid.cells_tot[i].state_consv_center[state_var] + grid.cells_tot[i].rhs_state[state_var]*dt;
            }
        }

        writeln("After update by ", dt);
        write("[");
        for i in grid.indicesInner {
            if !grid.cells_tot[i].solve_flag then continue;
            for state_var in grid.cells_tot[i].states_count {
                write(grid.cells_tot[i].state_consv_center[state_var], ", ");
            }
        }
        writeln("\b\b]");
    }

    proc updateCells (grid: borrowed Grid(?), dt: real(64)): void {
        if time_integration=="euler" then
                timeStepEuler(grid, dt);
        else if time_integration=="rk2" then
                timeStepRK2(grid, dt);
        else {
            writeln("Problem: ", time_integration, " time integration is un-supported!");
                exit(1);
        }
        consv_to_prims (grid);
        sync grid.cells_tot.updateFluff();
    }
}