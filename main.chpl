use Help;
use structures;
use reconstruction;
use conversion;
use initialize;
use boundary;

config const npoints: int(64)  = 5;
config const xmin:    real(64) = 0.0;
config const xmax:    real(64) = 5.0;
config const reconstruction_type: string = "constant";
config const cfl: real(64) = 0.3;
config const dt_ini: real(64) = 1.0e-04;
config const t_start: real(64) = 0.0;
config const t_stop: real(64) = 1.0;

proc main(args: [] string) {  
    var w: Wall;
    const nghosts: int(64) = compute_nghost(reconstruction_type);
    var grid = new owned Grid(xmin, xmax, npoints, nghosts);
    var computeDomain: domain(1) = {grid.indicesInner.low..grid.indicesInner.high};

    if debug { 
      writeln("cell centers: ");
      for i in grid.cells_tot.domain {
        write(grid.cells_tot[i].center);
        write(" ");
      }
      writeln();
      writeln("left wall positions: ");
      for i in grid.cells_tot.domain {
        write(grid.cells_tot[i].wall_left.position);
        write(" ");
      }
      writeln();
      writeln("right wall positions: ");
      for i in grid.cells_tot.domain {
        write(grid.cells_tot[i].wall_right.position);
        write(" ");
      }
      writeln();
    }

    var stepNumber: int(64) = 0;
    var time: real(64) = t_start;
    var delta_t: real(64) = dt_ini;
    var cell_size_min: real(64) = min reduce [c in grid.cells_tot] c.cell_size;
    if debug then writeln("Min cell size: ", cell_size_min);

    init_field(grid);
    prims_to_consv(grid);
    set_boundary(grid);
    interpolate_edges(grid);

    while time<t_stop do {

      stepNumber += 1;
      time += delta_t;
    }

    writeln("Cells: ");
    for i in grid.indicesAll {
      var tag: string = "";
      if !computeDomain.contains(i) then tag = ": boundary";
      writeln("x = ", grid.cells_tot[i].center, ": (", grid.cells_tot[i].state_consv_center, ", ", grid.cells_tot[i].state_prims_center, ")", tag);
    }
    writeln("Walls: ");
    for i in grid.indicesAllStag {
      var tag: string = "";
      if i==grid.indicesAllStag.low || i==grid.indicesAllStag.high then tag = ": boundary";
      writeln("x = ", grid.walls_tot[i].position, ": ", grid.walls_tot[i].state_consv_left, " | ", grid.walls_tot[i].state_consv_right, tag);
  }

    /*
    writeln("Left boundary: ");
    forall (i, j) in grid.indicesAll.boundaries(0, -1) {
      writeln("element ", i, " side ", j);
    }
    writeln("Right boundary: ");
    forall (i, j) in grid.indicesAll.boundaries(0, 1) {
      writeln("element ", i, " side ", j);
    }
    */
}
