use Help;
use Math;
use utilities;
use structures;
use reconstruction;
use conversion;
use initialize;
use boundary;
use riemann;
use time_stepper;

config const npoints: int(64)  = 5;
config const xmin:    real(64) = 0.0;
config const xmax:    real(64) = 5.0;
config const reconstruction_type: string = "constant";
config const cfl: real(64) = 0.3;
config const dt_ini: real(64) = 1.0e-04;
config const t_start: real(64) = 0.0;
config const t_stop: real(64) = 3.0;
config const vel_adv: real(64) = 1.0;
config const dt_max: real(64) = 1.0e-02;
config const max_steps: int(64) = 10000000000000;
config const freq: int(64) = 20;
config const output_interval: real(64) = 0.1;

proc shouldOutput(time: real(64), output_dt: real(64), dt: real(64)): bool {
  return abs(mod(time, output_dt)) < 1.1*dt;
}

proc updateCells (grid: borrowed Grid(?), dt: real(64)): void {
  var computeDomain: domain(1) = {grid.indicesInner.low..grid.indicesInner.high};
  solve_at_walls(grid, vel_adv);
  sync forall i in grid.indicesInner {
    if !grid.cells_tot[i].solve_flag then continue;
    var rhs_state: [grid.cells_tot[i].states_count] real(64);
    var consv_updated: [grid.cells_tot[i].states_count] real(64);
    for state_var in grid.cells_tot[i].states_count {
      rhs_state[state_var] = (grid.cells_tot[i].wall_left.flux_solve[state_var] - grid.cells_tot[i].wall_right.flux_solve[state_var])/grid.cells_tot[i].cell_size;
    }
    consv_updated = timeStep(grid.cells_tot[i].state_consv_center, rhs_state, dt);
    for state_var in grid.cells_tot[i].states_count do 
      grid.cells_tot[i].state_consv_center[state_var] = consv_updated[state_var];
  }
  consv_to_prims (grid);
  sync grid.cells_tot.updateFluff();
}

proc main(args: [] string) {  
    var w: Wall;
    const nghosts: int(64) = compute_nghost(reconstruction_type);
    var grid = new owned Grid(xmin, xmax, npoints, nghosts);
    var computeDomain: domain(1) = {grid.indicesInner.low..grid.indicesInner.high};

    /*
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
    }*/

    var stepNumber: int(64) = 0;
    var time: real(64) = t_start;
    var delta_t: real(64) = dt_ini;
    var cell_size_min: real(64) = min reduce [c in grid.cells_tot] c.cell_size;
    if debug then writeln("Min cell size: ", cell_size_min);

    init_field(grid);
    prims_to_consv(grid);
    set_boundary(grid);
    interpolate_edges(grid);

    /*
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
    */
    var positions: [computeDomain] real(64);
    var state: [computeDomain] real(64);
    sync forall i in computeDomain {
      positions[i] = grid.cells_tot[i].center;
      state[i] = grid.cells_tot[i].state_prims_center[grid.cells_tot[i].states_count.low];
    }
    writeArraysToFile("./output."+stepNumber:string+".txt", positions, state);

    if delta_t>dt_max then delta_t = dt_max;
    writeln("Starting computation ...");

    writeln("time ", time, " (step ", stepNumber ,"): dt = ", delta_t);
    while time<t_stop && stepNumber<max_steps do {
      /* computation for each loop starts here */
      updateCells(grid, delta_t);
      // XXX: Check if this is needed later 
      // prims_to_consv(grid);
      set_boundary(grid);
      interpolate_edges(grid);
      /* computation for each loop ends here */
      stepNumber += 1;
      time += delta_t;
      delta_t = cfl*cell_size_min/vel_adv;
      if delta_t>dt_max then delta_t = dt_max;
      if stepNumber%freq==0 || stepNumber==max_steps then 
        writeln("time ", time, " (step ", stepNumber ,"): dt = ", delta_t);
      if shouldOutput(time, output_interval, delta_t) && (stepNumber>1 || max_steps<10000000000000) {
        writeln("Dumping ./output.", stepNumber, ".txt to disk");
        sync forall i in computeDomain {
          positions[i] = grid.cells_tot[i].center;
          state[i] = grid.cells_tot[i].state_prims_center[grid.cells_tot[i].states_count.low];
        }
        writeArraysToFile("./output."+stepNumber:string+".txt", positions, state);
      }
    }
    sync forall i in computeDomain {
      positions[i] = grid.cells_tot[i].center;
      state[i] = grid.cells_tot[i].state_prims_center[grid.cells_tot[i].states_count.low];
    }
    writeArraysToFile("./output."+stepNumber:string+".txt", positions, state);

    /*
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
    */
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
