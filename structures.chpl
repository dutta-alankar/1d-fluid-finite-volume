module structures {
    import StencilDist.stencilDist;
    import utilities;

    config const debug: bool = false;
    config const max_states: int(64) = 5;
    const state_domain: domain(1) = {1..max_states};

    record Wall {
        var position: real(64);
        var states_count: domain(1) = state_domain;
        var state_consv_left: [state_domain] real(64);
        var state_consv_right: [state_domain] real(64);
        var state_consv_solve: [state_domain] real(64);
        var flux_solve: [state_domain] real(64);
        var dummy_initialized: bool;

        proc init (): void {
            this.position = 0.0;
            this.dummy_initialized = true;
            init this;
            for state_var in this.states_count {
                this.state_consv_left[state_var]  = -state_var: real(64);
                this.state_consv_right[state_var] = -state_var: real(64);
                this.state_consv_solve[state_var] = -state_var: real(64);
                this.flux_solve[state_var] = -state_var: real(64);
            }
        }

        proc init (position: real(64)): void {
            init();
            this.position = position;
            this.dummy_initialized = false;
        }

        // XXX: add reconstruction to fill state
    }

    record Cell {
        var wall_left: Wall;
        var wall_right: Wall;
        var indices: int(64);
        var center: real(64);
        var cell_size: real(64);
        var states_count: domain(1) = state_domain;
        var state_consv_center: [state_domain] real(64);
        var state_prims_center: [state_domain] real(64);
        var dummy_initialized: bool;

        proc init (): void {
            this.wall_left  = new Wall();
            this.wall_right = new Wall();
            this.indices = 0;
            this.center = 0.0;
            this.cell_size = 0.0;
            this.dummy_initialized = true;
            init this;
            for state_var in this.states_count {
                this.state_consv_center[state_var] = -state_var: real(64);
                this.state_prims_center[state_var] = -state_var: real(64);
            }
        }

        proc init (const ref wall_left: Wall, const ref wall_right: Wall, indices: int(64)): void {
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

    class Grid {
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
        var indicesInner, indicesAll, indicesAllStag: domain(?);
        var cells_tot: [indicesAll]  Cell;
        var walls_tot: [indicesAllStag] Wall;

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
            var indicesInnerDomain: domain(1) = {1..this.npoints_int};
            this.indicesInner = indicesInnerDomain dmapped new stencilDist(indicesInnerDomain, fluff=(this.nghosts,), periodic=false);
            this.indicesAll = indicesInnerDomain.expand((this.nghosts,));
            var indicesAllStagDomain: domain(1) = {this.indicesAll.low..this.indicesAll.high+1};
            this.indicesAllStag = indicesAllStagDomain dmapped new stencilDist(indicesAllStagDomain, fluff=(this.nghosts,), periodic=false);
            init this; 
            this.create_grid(); 
        }

        proc create_grid (): void {
            // uniform grid
            var dx: real(64) = (this.xmax-this.xmin)/this.npoints_int;
            var x_left:  [this.indicesAll] real(64) = utilities.linspace(this.xmin-this.nghosts*dx, this.xmax+(this.nghosts-1)*dx, this.npoints_tot, this.indicesAll);
            var x_right: [this.indicesAll] real(64) = x_left + dx;
            if debug {
                writeln("indicesAll: ", this.indicesAll);
                writeln("indicesInner: ", this.indicesInner);
                writeln("indicesAllStag: ", this.indicesAllStag);
                writeln("Left wall position info for debugging: ");
                writeln(x_left);
                writeln("Right wall position info for debugging: ");
                writeln(x_right);
            }
            // create the wall
            sync {
                forall i in this.indicesAll do
                    this.walls_tot[i] = new Wall(x_left[i]);
                this.walls_tot[this.indicesAllStag.high] = new Wall(x_right[this.indicesAll.high]);
            }
            sync this.walls_tot.updateFluff();
            if debug {
                writeln("Wall info for debugging: ");
                for i in this.walls_tot.domain {
                    if !this.walls_tot[i].dummy_initialized then
                        writeln("i=", i, " -> pos=", this.walls_tot[i].position, " at locale ", this.walls_tot[i].locale.id);
                    else
                        writeln("i=", i, " -> dummy!", " at locale ", this.walls_tot[i].locale.id);
                }
            }
            // create the cells
            sync forall i in this.indicesAll {
                assert(!this.walls_tot[i].dummy_initialized, "Problem: Wall "+i:string+" is dummy");
                assert(!this.walls_tot[i+1].dummy_initialized, "Problem: Wall "+(i+1):string+" is dummy");
                this.cells_tot[i] = new Cell(this.walls_tot[i], this.walls_tot[i+1], i); 
            }
            sync this.cells_tot.updateFluff();
        }
    }
}