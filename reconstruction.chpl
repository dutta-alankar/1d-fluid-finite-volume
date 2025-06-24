module reconstruction {
    use structures;

    var reconstruction_coeff_domain: domain(1);
    var reconstruction_coeff: [reconstruction_coeff_domain] real(64);
    
    proc compute_nghost(const ref reconstruction_type: string ): int(64) {
        if reconstruction_type=="constant" {
            reconstruction_coeff_domain = {1..1};
            for i in reconstruction_coeff_domain {
                reconstruction_coeff[i] = 1.0;
            }
            return 1;
        }
        else if reconstruction_type=="linear" {
            reconstruction_coeff_domain = {1..2};
            for i in reconstruction_coeff_domain {
                reconstruction_coeff[i] = 1.0;
            }
            return 2;
        }
        else {
            writeln("Problem: Invaid reconstruction scheme!");
            reconstruction_coeff_domain = {1..1};
            for i in reconstruction_coeff_domain {
                reconstruction_coeff[i] = 1.0;
            }
            return 1;
        }
    }

    proc interpolate_edges (grid: borrowed Grid(?)): void {
        // TODO: might need indicesInner
        sync forall i in grid.indicesAll {
            var sgn_dist_left: real(64)  = -(grid.cells_tot[i].center-grid.cells_tot[i].wall_left.position)/grid.cells_tot[i].cell_size;
            var sgn_dist_right: real(64) = -(grid.cells_tot[i].center-grid.cells_tot[i].wall_right.position)/grid.cells_tot[i].cell_size;
            var pow: int(64) = 0;
            for j in reconstruction_coeff.domain {
                for state_var in grid.cells_tot[i].state_consv_center.domain {
                    grid.cells_tot[i].wall_left.state_consv_right =  grid.cells_tot[i].state_consv_center[state_var] * reconstruction_coeff[j] * sgn_dist_left**pow;
                    grid.cells_tot[i].wall_right.state_consv_left =  grid.cells_tot[i].state_consv_center[state_var] * reconstruction_coeff[j] * sgn_dist_right**pow;
                }
                pow = pow+1;
            }
        }
        sync grid.walls_tot.updateFluff();
    }
}