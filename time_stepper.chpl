module time_stepper {
    use globals;
    use structures;
    use riemann;
    use conversion;

    proc timeStepEuler (grid: borrowed Grid(?), dt: real(64)): void {
        /* explicit Euler first order time integration */
        // riemann problem solve
        solve_at_walls(grid, (vel_adv_x, vel_adv_y, vel_adv_z));
        var computeDomain: domain(3) = {grid.indicesInner.low[0]..grid.indicesInner.high[0], grid.indicesInner.low[1]..grid.indicesInner.high[1], grid.indicesInner.low[2]..grid.indicesInner.high[2]};
        forall (i,j,k) in grid.indicesInner {
            if !grid.cells_tot[i,j,k].solve_flag then {
                continue;
            }
            for state_var in grid.cells_tot[i,j,k].states_count {
                grid.cells_tot[i,j,k].rhs_state[state_var] = (
                    (grid.cells_tot[i,j,k].wall_left.flux_solve[state_var] - grid.cells_tot[i,j,k].wall_right.flux_solve[state_var])/grid.cells_tot[i,j,k].cell_size[0] +
                    (grid.cells_tot[i,j,k].wall_front.flux_solve[state_var] - grid.cells_tot[i,j,k].wall_back.flux_solve[state_var])/grid.cells_tot[i,j,k].cell_size[1] +
                    (grid.cells_tot[i,j,k].wall_bottom.flux_solve[state_var] - grid.cells_tot[i,j,k].wall_top.flux_solve[state_var])/grid.cells_tot[i,j,k].cell_size[2]
                );
            }
        }
        forall (i,j,k) in grid.indicesInner {
            if !grid.cells_tot[i,j,k].solve_flag then {
                continue;
            }
            for state_var in grid.cells_tot[i,j,k].states_count {
                grid.cells_tot[i,j,k].state_consv_center[state_var] = grid.cells_tot[i,j,k].state_consv_center[state_var] + grid.cells_tot[i,j,k].rhs_state[state_var]*dt;
            }
        }
    }

    proc timeStepRK2 (grid: borrowed Grid(?), dt: real(64)): void {
        /* explicit RK2 second order time integration */
        var grid_tmp: Grid(?) = new Grid();
        // RK2 step 1 : use the already calculated flux for half timeStep
        grid.deepCopy(grid_tmp);
        timeStepEuler(grid_tmp, 0.5*dt);
        consv_to_prims(grid_tmp);
        // RK2 step 2 : use the half timeStep conservative states to reconstruct and evaluate full time state
        set_boundary(grid_tmp);
        interpolate_edges(grid_tmp);
        // riemann problem solve
        solve_at_walls(grid_tmp, (vel_adv_x, vel_adv_y, vel_adv_z));
        // set the rhs according to flux at half-time step
        // XXX: state_consv_solve & flux_solve not copied from grid_tmp (might be needed later)
        var computeDomain: domain(3) = {grid.indicesInner.low[0]..grid.indicesInner.high[0], grid.indicesInner.low[1]..grid.indicesInner.high[1], grid.indicesInner.low[2]..grid.indicesInner.high[2]};
        forall (i,j,k) in grid.indicesInner {
            if !grid.cells_tot[i,j,k].solve_flag then {
                continue;
            }
            for state_var in grid.cells_tot[i,j,k].states_count {
                grid.cells_tot[i,j,k].rhs_state[state_var] = (
                    (grid_tmp.cells_tot[i,j,k].wall_left.flux_solve[state_var] - grid_tmp.cells_tot[i,j,k].wall_right.flux_solve[state_var])/grid_tmp.cells_tot[i,j,k].cell_size[0] +
                    (grid_tmp.cells_tot[i,j,k].wall_front.flux_solve[state_var] - grid_tmp.cells_tot[i,j,k].wall_back.flux_solve[state_var])/grid_tmp.cells_tot[i,j,k].cell_size[1] +
                    (grid_tmp.cells_tot[i,j,k].wall_bottom.flux_solve[state_var] - grid_tmp.cells_tot[i,j,k].wall_top.flux_solve[state_var])/grid_tmp.cells_tot[i,j,k].cell_size[2]
                );
            }
        }
        // use the rhs fluxes based on half timestep to evolve conservative states
        forall (i,j,k) in grid.indicesInner {
            if !grid.cells_tot[i,j,k].solve_flag then {
                continue;
            }
            for state_var in grid.cells_tot[i,j,k].states_count {
                grid.cells_tot[i,j,k].state_consv_center[state_var] = grid.cells_tot[i,j,k].state_consv_center[state_var] + grid.cells_tot[i,j,k].rhs_state[state_var]*dt;
            }
        }
    }

    proc updateCells (grid: borrowed Grid(?), dt: real(64)): void {
        if time_integration=="euler" then {
                timeStepEuler(grid, dt);
        } else if time_integration=="rk2" then {
                timeStepRK2(grid, dt);
        } else {
            writeln("Problem: ", time_integration, " time integration is un-supported!");
            exit(1);
        }
        consv_to_prims (grid);
        sync grid.cells_tot.updateFluff();
    }
}