module structures {
    import StencilDist.stencilDist;
    import utilities;

    const state_domain: domain(1) = {1..5};

    record wall {
        var position: real(64);
        var states_count: domain(1);
        var state_prims: [states_count] real(64);
        var state_consv: [states_count] real(64);
        var dummy_initialized: bool;

        proc init (): void {
            this.position = 0.0;
            this.states_count = state_domain;
            this.state_prims = [i in this.states_count] 0.0;
            this.state_consv = [i in this.states_count] 0.0;
            this.dummy_initialized = true;
        }

        proc init (position: real(64)): void {
            init();
            this.position = position;
            this.dummy_initialized = false;
        }

        // XXX: add reconstruction to fill state
    }

    record cell {
        var wall_left: wall;
        var wall_right: wall;
        var indices: int(64);
        var center: real(64);
        var cell_size: real(64);
        var dummy_initialized: bool;

        proc init (): void {
            this.wall_left  = new wall();
            this.wall_right = new wall();
            this.indices = 0;
            this.center = 0.0;
            this.cell_size = 0.0;
            this.dummy_initialized = true;
        }

        proc init (wall_left: wall, wall_right: wall, indices: int(64)): void {
            this.wall_left  = wall_left;
            this.wall_right = wall_right;
            this.indices = indices;
            var xleft: real(64)  = this.wall_left.position;
            var xright: real(64) = this.wall_right.position;
            this.center = 0.5*(xleft+xright);
            this.cell_size = xright-xleft;
            this.dummy_initialized = false;
        }
    }

    class grid {
        /* In my code all indices start from 1 */
        var xmin: real(64);
        var xmax: real(64);
        var nghosts: int(64);
        var npoints_int: int(64);
        var indx_beg_int: int(64);
        var indx_end_int: int(64);
        var indx_beg_tot: int(64);
        var indx_end_tot: int(64);
        var npoints_tot: int(64);
        // initialize with a dummy
        var indicesAll, indicesInner, indicesStag: domain(?);
        var cells_tot: [indicesAll]  cell;
        var walls_tot: [indicesStag] wall;

        proc init (xmin: real(64), xmax: real(64), npoints: int(64), nghosts: int(64)): void {
            this.xmin = xmin;
            this.xmax = xmax;
            this.nghosts = nghosts; 
            this.npoints_int = npoints;           
            this.indx_beg_int = nghosts+1;
            this.indx_end_int = npoints_int + this.indx_beg_int;
            this.indx_beg_tot = 1;
            this.indx_end_tot = npoints_int + 2*nghosts;
            this.npoints_tot  = this.indx_end_tot;
            this.indicesAll = stencilDist.createDomain({1..this.npoints_tot}, fluff=(this.nghosts,));
            this.indicesInner = this.indicesAll.expand((-this.nghosts,));
            this.indicesStag = stencilDist.createDomain({1..this.npoints_tot+1}, fluff=(this.nghosts,));
            init this; 
            this.create_grid(); 
        }

        proc create_grid (): void {
            // uniform grid
            var dx: real(64) = (this.xmax-this.xmin)/this.npoints_int;
            var x_left:  [this.indicesAll] real(64) = utilities.linspace(this.xmin-this.nghosts*dx, this.xmax+(this.nghosts-1)*dx, this.npoints_tot, this.indicesAll);
            var x_right: [this.indicesAll] real(64) = x_left + dx;
            // create the wall
            sync {
                forall i in this.indicesAll {
                    this.walls_tot[i] = new wall(x_left[i]);
                }
                this.walls_tot[this.indicesStag.high] = new wall(x_right[this.indicesAll.high]);
            }
            this.walls_tot.updateFluff();
            writeln("Walls on locale ", here.id, ":");
            forall i in this.walls_tot.domain {
                if !this.walls_tot[i].dummy_initialized then
                    writeln("i=", i, " -> pos=", this.walls_tot[i].position, " at locale ", this.walls_tot[i].locale.id);
                else
                    writeln("i=", i, " -> dummy!", " at locale ", this.walls_tot[i].locale.id);
            }
            // create the cells
            sync {
                forall i in this.indicesAll {
                    // assert(!this.walls_tot[i].dummy_initialized, "Problem: Wall "+i:string+" is dummy");
                    // assert(!this.walls_tot[i+1].dummy_initialized, "Problem: Wall "+(i+1):string+" is dummy");
                    this.cells_tot[i] = new cell(this.walls_tot[i], this.walls_tot[i+1], i); 
                }
            }
            this.cells_tot.updateFluff();

        }
    }
}