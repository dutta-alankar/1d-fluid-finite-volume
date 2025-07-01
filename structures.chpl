module structures {
    use globals;
    import StencilDist.stencilDist;
    import utilities;

    config const debug: bool = false;
    config const max_states: int(64) = 5;
    const state_domain = {1..max_states};

    record Wall {
        var position: (real(64), real(64), real(64));
        var states_count: domain(1) = state_domain;
        var state_consv_left: [state_domain] real(64);
        var state_consv_right: [state_domain] real(64);
        var state_consv_front: [state_domain] real(64);
        var state_consv_back: [state_domain] real(64);
        var state_consv_bottom: [state_domain] real(64);
        var state_consv_top: [state_domain] real(64);
        var state_consv_solve: [state_domain] real(64);
        var flux_solve: [state_domain] real(64);
        var dummy_initialized: bool;

        proc init (): void {
            this.position = (0.0, 0.0, 0.0);
            this.dummy_initialized = true;
            init this;
            for state_var in this.states_count {
                this.state_consv_left[state_var]  = -state_var: real(64);
                this.state_consv_right[state_var] = -state_var: real(64);
                this.state_consv_front[state_var] = -state_var: real(64);
                this.state_consv_back[state_var] = -state_var: real(64);
                this.state_consv_bottom[state_var] = -state_var: real(64);
                this.state_consv_top[state_var] = -state_var: real(64);
                this.state_consv_solve[state_var] = -state_var: real(64);
                this.flux_solve[state_var] = -state_var: real(64);
            }
        }

        proc init (position: (real(64), real(64), real(64))): void {
            init();
            this.position = position;
            this.dummy_initialized = false;
        }

        // XXX: add reconstruction to fill state
    }

    record Cell {
        var wall_left: Wall;
        var wall_right: Wall;
        var wall_front: Wall;
        var wall_back: Wall;
        var wall_bottom: Wall;
        var wall_top: Wall;
        var indices: (int(64), int(64), int(64));
        var center: (real(64), real(64), real(64));
        var cell_size: (real(64), real(64), real(64));
        var states_count: domain(1) = state_domain;
        var state_consv_center: [state_domain] real(64);
        var state_prims_center: [state_domain] real(64);
        var rhs_state: [state_domain] real(64);
        var dummy_initialized: bool;
        var solve_flag: bool;

        proc init (): void {
            this.wall_left  = new Wall();
            this.wall_right = new Wall();
            this.wall_front = new Wall();
            this.wall_back = new Wall();
            this.wall_bottom = new Wall();
            this.wall_top = new Wall();
            this.indices = (0, 0, 0);
            this.center = (0.0, 0.0, 0.0);
            this.cell_size = (0.0, 0.0, 0.0);
            this.dummy_initialized = true;
            this.solve_flag = true;
            init this;
            for state_var in this.states_count {
                this.state_consv_center[state_var] = -state_var: real(64);
                this.state_prims_center[state_var] = -state_var: real(64);
                this.rhs_state[state_var] = -state_var: real(64);
            }
        }

        proc init (const ref wall_left: Wall, const ref wall_right: Wall,
                   const ref wall_front: Wall, const ref wall_back: Wall,
                   const ref wall_bottom: Wall, const ref wall_top: Wall,
                   indices: (int(64), int(64), int(64))): void {
            this.wall_left  = wall_left;
            this.wall_right = wall_right;
            this.wall_front = wall_front;
            this.wall_back = wall_back;
            this.wall_bottom = wall_bottom;
            this.wall_top = wall_top;
            this.indices = indices;
            var xleft: real(64)  = this.wall_left.position[0];
            var xright: real(64) = this.wall_right.position[0];
            var yfront: real(64) = this.wall_front.position[1];
            var yback: real(64)  = this.wall_back.position[1];
            var zbottom: real(64) = this.wall_bottom.position[2];
            var ztop: real(64)   = this.wall_top.position[2];
            this.center = (0.5*(xleft+xright), 0.5*(yfront+yback), 0.5*(zbottom+ztop));
            this.cell_size = (xright-xleft, yback-yfront, ztop-zbottom);
            this.dummy_initialized = false;
            this.solve_flag = true;
        }
    }

    class Grid {
        /* In my code all indices start from 1 */
        var xmin: (real(64), real(64), real(64));
        var xmax: (real(64), real(64), real(64));
        var nghosts: (int(64), int(64), int(64));
        var npoints_int: (int(64), int(64), int(64));
        var indx_beg_int: (int(64), int(64), int(64));
        var indx_end_int: (int(64), int(64), int(64));
        var indx_beg_tot: (int(64), int(64), int(64));
        var indx_end_tot: (int(64), int(64), int(64));
        var npoints_tot: (int(64), int(64), int(64));
        // initialize with a dummy
        var indicesInner, indicesAll, indicesAllStag: domain(?);
        var cells_tot: [indicesAll]  Cell;
        var walls_tot: [indicesAllStag] Wall;

        proc init (xmin: (real(64), real(64), real(64)) = (0.0, 0.0, 0.0), xmax: (real(64), real(64), real(64)) = (1.0, 1.0, 1.0), npoints: (int(64), int(64), int(64)) = (1, 1, 1), nghosts: (int(64), int(64), int(64)) = (1, 1, 1)): void {
            this.xmin = xmin;
            this.xmax = xmax;
            this.nghosts = nghosts; 
            this.npoints_int = npoints;           
            this.indx_beg_int = (nghosts[0]+1, nghosts[1]+1, nghosts[2]+1);
            this.indx_end_int = (npoints_int[0] + this.indx_beg_int[0], npoints_int[1] + this.indx_beg_int[1], npoints_int[2] + this.indx_beg_int[2]);
            this.indx_beg_tot = (1, 1, 1);
            this.indx_end_tot = (npoints_int[0] + 2*nghosts[0], npoints_int[1] + 2*nghosts[1], npoints_int[2] + 2*nghosts[2]);
            this.npoints_tot  = this.indx_end_tot;
            var indicesInnerDomain: domain(3) = {1..this.npoints_int[0], 1..this.npoints_int[1], 1..this.npoints_int[2]};
            this.indicesInner = indicesInnerDomain dmapped new stencilDist(indicesInnerDomain, fluff=(this.nghosts[0], this.nghosts[1], this.nghosts[2]), periodic=false);
            this.indicesAll = indicesInnerDomain.expand((this.nghosts[0], this.nghosts[1], this.nghosts[2]));
            var indicesAllStagDomain: domain(3) = {this.indicesAll.low[0]..this.indicesAll.high[0]+1, this.indicesAll.low[1]..this.indicesAll.high[1]+1, this.indicesAll.low[2]..this.indicesAll.high[2]+1};
            this.indicesAllStag = indicesAllStagDomain dmapped new stencilDist(indicesAllStagDomain, fluff=(this.nghosts[0], this.nghosts[1], this.nghosts[2]), periodic=false);
            init this;
            if npoints[0]>1 || npoints[1]>1 || npoints[2]>1 then
                this.create_grid(); 
        }

        proc deepCopy (grid: borrowed Grid(?)): void {
            grid.xmin = this.xmin;
            grid.xmax = this.xmax;
            grid.nghosts = this.nghosts; 
            grid.npoints_int = this.npoints_int;           
            grid.indx_beg_int = (this.nghosts[0]+1, this.nghosts[1]+1, this.nghosts[2]+1);
            grid.indx_end_int = (this.npoints_int[0] + this.indx_beg_int[0], this.npoints_int[1] + this.indx_beg_int[1], this.npoints_int[2] + this.indx_beg_int[2]);
            grid.indx_beg_tot = (1, 1, 1);
            grid.indx_end_tot = (this.npoints_int[0] + 2*this.nghosts[0], this.npoints_int[1] + 2*this.nghosts[1], this.npoints_int[2] + 2*this.nghosts[2]);
            grid.npoints_tot  = grid.indx_end_tot;
            var indicesInnerDomain: domain(3) = {1..this.npoints_int[0], 1..this.npoints_int[1], 1..this.npoints_int[2]};
            grid.indicesInner = indicesInnerDomain dmapped new stencilDist(indicesInnerDomain, fluff=(this.nghosts[0], this.nghosts[1], this.nghosts[2]), periodic=false);
            grid.indicesAll = indicesInnerDomain.expand((this.nghosts[0], this.nghosts[1], this.nghosts[2]));
            var indicesAllStagDomain: domain(3) = {this.indicesAll.low[0]..this.indicesAll.high[0]+1, this.indicesAll.low[1]..this.indicesAll.high[1]+1, this.indicesAll.low[2]..this.indicesAll.high[2]+1};
            grid.indicesAllStag = indicesAllStagDomain dmapped new stencilDist(indicesAllStagDomain, fluff=(this.nghosts[0], this.nghosts[1], this.nghosts[2]), periodic=false);
            forall i in grid.indicesAllStag do
                grid.walls_tot[i] = this.walls_tot[i];
            forall i in grid.indicesAll do
                grid.cells_tot[i] = this.cells_tot[i];
            // XXX: might not need this and removing this can improve performance
            sync grid.walls_tot.updateFluff();
            sync grid.cells_tot.updateFluff();
        }

        proc create_grid (): void {
            // uniform grid
            var dx: real(64) = (this.xmax[0]-this.xmin[0])/this.npoints_int[0];
            var dy: real(64) = (this.xmax[1]-this.xmin[1])/this.npoints_int[1];
            var dz: real(64) = (this.xmax[2]-this.xmin[2])/this.npoints_int[2];

            var x_left:  [this.indicesAll] real(64);
            var x_right: [this.indicesAll] real(64);
            var y_front: [this.indicesAll] real(64);
            var y_back:  [this.indicesAll] real(64);
            var z_bottom: [this.indicesAll] real(64);
            var z_top:   [this.indicesAll] real(64);

            forall (i,j,k) in this.indicesAll {
                x_left[i,j,k] = this.xmin[0]-this.nghosts[0]*dx + (i-1)*dx;
                x_right[i,j,k] = x_left[i,j,k] + dx;
                y_front[i,j,k] = this.xmin[1]-this.nghosts[1]*dy + (j-1)*dy;
                y_back[i,j,k] = y_front[i,j,k] + dy;
                z_bottom[i,j,k] = this.xmin[2]-this.nghosts[2]*dz + (k-1)*dz;
                z_top[i,j,k] = z_bottom[i,j,k] + dz;
            }

            var computeDomain: domain(3) = {this.indicesInner.low[0]..this.indicesInner.high[0], this.indicesInner.low[1]..this.indicesInner.high[1], this.indicesInner.low[2]..this.indicesInner.high[2]};

            // create the walls
            sync {
                forall (i,j,k) in this.indicesAll {
                    this.walls_tot[i,j,k] = new Wall((x_left[i,j,k], y_front[i,j,k], z_bottom[i,j,k]));
                }
                this.walls_tot[this.indicesAllStag.high[0], this.indicesAllStag.high[1], this.indicesAllStag.high[2]] = new Wall((x_right[this.indicesAll.high[0], this.indicesAll.high[1], this.indicesAll.high[2]], y_back[this.indicesAll.high[0], this.indicesAll.high[1], this.indicesAll.high[2]], z_top[this.indicesAll.high[0], this.indicesAll.high[1], this.indicesAll.high[2]]));
            }
            sync this.walls_tot.updateFluff();

            // create the cells
            forall (i,j,k) in this.indicesAll {
                assert(!this.walls_tot[i,j,k].dummy_initialized, "Problem: Wall "+(i,j,k):string+" is dummy");
                assert(!this.walls_tot[i+1,j,k].dummy_initialized, "Problem: Wall "+(i+1,j,k):string+" is dummy");
                assert(!this.walls_tot[i,j+1,k].dummy_initialized, "Problem: Wall "+(i,j+1,k):string+" is dummy");
                assert(!this.walls_tot[i,j,k+1].dummy_initialized, "Problem: Wall "+(i,j,k+1):string+" is dummy");
                this.cells_tot[i,j,k] = new Cell(this.walls_tot[i,j,k], this.walls_tot[i+1,j,k],
                                                this.walls_tot[i,j,k], this.walls_tot[i,j+1,k],
                                                this.walls_tot[i,j,k], this.walls_tot[i,j,k+1],
                                                (i,j,k));
                if !computeDomain.contains((i,j,k)) then this.cells_tot[i,j,k].solve_flag = false;
            }
            sync this.cells_tot.updateFluff();
        }
    }
}