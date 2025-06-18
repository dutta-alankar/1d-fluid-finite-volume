module structures {
    import StencilDist.stencilDist;
    import utilities;

    const state_domain: domain(1) = {1..5};

    class wall {
        var position: real(64);
        var state: [state_domain] real(64);

        proc init (position: real(64)): void {
            this.position = position;
        }

        // XXX: add reconstruction to fill state
    }

    class cell {
        var left_wall: owned wall;
        var right_wall: owned wall;
        var indices: int(64);
        var center: real(64);
        var cell_size: real(64);

        proc init (xleft: real(64), xright: real(64), indices: int(64)): void {
            this.left_wall  = new owned wall(xleft);
            this.right_wall = new owned wall(xright);
            this.indices = indices;
            this.center = 0.5*(this.left_wall.position+this.right_wall.position);
            this.cell_size = this.right_wall.position-this.left_wall.position;
        }
    }
    /*
    class wall {
        const cell_left: borrowed cell;
        const cell_right: borrowed cell;

        proc init (cell_left: borrowed cell, cell_right: borrowed cell) {
            this.cell_left = cell_left;
            this.cell_right = cell_right;
        }
    }
    */
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
        var indices: domain(?);
        var indicesInner: domain(?);
        var indicesStag: domain(?);
        var cells_tot: [indices] owned cell?;

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
            this.indices = stencilDist.createDomain({1..this.npoints_tot}, fluff=(this.nghosts,));
            this.indicesInner = indices.expand((-this.nghosts,));
            this.indicesStag = stencilDist.createDomain({1..this.npoints_tot+1}, fluff=(this.nghosts,));
            init this; 
            this.create_grid(); 
        }

        proc create_grid (): void {
            // uniform grid
            var dx: real(64) = (this.xmax-this.xmin)/this.npoints_int;
            var x_left:  [this.indices] real(64) = utilities.linspace(this.xmin-this.nghosts*dx, this.xmax+(this.nghosts-1)*dx, this.npoints_tot, this.indices);
            var x_right: [this.indices] real(64) = x_left + dx;
            forall i in this.indices {
                this.cells_tot[i] = new owned cell(x_left[i], x_right[i], i);
            }
        }

    }
}