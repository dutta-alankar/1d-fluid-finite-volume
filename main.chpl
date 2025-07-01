use Help;
use utilities;
use structures;
use reconstruction;
use conversion;
use initialize;
use boundary;
use riemann;
use time_stepper;
use dump;
use logo;
use globals;

proc warnings (): void {
  if time_integration=="euler" && reconstruction_type=="linear" then
    writeln("WARNING: FTCS scheme is unconditionally unstable.");
}

proc prepare_run (grid: borrowed Grid(?)): void {
  sync {
    init_field(grid);
    set_boundary(grid);
    prims_to_consv(grid);
  }
  interpolate_edges(grid);
  // solve_at_walls(grid, vel_adv);
}

proc main(args: [] string) { 
    print_logo();
    warnings();
    var w: Wall;
    const nghosts: (int(64), int(64), int(64)) = compute_nghost(reconstruction_type);
    var grid = new owned Grid((xmin_x, xmin_y, xmin_z), (xmax_x, xmax_y, xmax_z), (npoints_x, npoints_y, npoints_z), nghosts);
    var computeDomain: domain(3) = {grid.indicesInner.low[0]..grid.indicesInner.high[0], grid.indicesInner.low[1]..grid.indicesInner.high[1], grid.indicesInner.low[2]..grid.indicesInner.high[2]};

    var stepNumber: int(64) = 0;
    var time: real(64) = t_start;
    var delta_t: real(64) = dt_ini;
    var cell_size_min: real(64) = min(min reduce [c in grid.cells_tot] c.cell_size[0], min reduce [c in grid.cells_tot] c.cell_size[1], min reduce [c in grid.cells_tot] c.cell_size[2]);
    if debug then writeln("Min cell size: ", cell_size_min, ", ghosts: ", nghosts);

    sync prepare_run(grid);
    dump_ascii_to_disk(grid, stepNumber, time);

    if delta_t>dt_max then delta_t = dt_max;
    if max_steps>0 then writeln("\n\nStarting computation ...");

    if max_steps>0 then writeln("time ", time, " (step ", stepNumber ,"): dt = ", delta_t);
    while time<t_stop && stepNumber<max_steps {
      /* computation for each loop starts here */
      updateCells(grid, delta_t);
      // XXX: Check if this is needed later 
      // prims_to_consv(grid);
      set_boundary(grid);
      interpolate_edges(grid);
      /* computation for each loop ends here */
      stepNumber += 1;
      time += delta_t;
      delta_t = cfl*cell_size_min/max(abs(vel_adv_x), abs(vel_adv_y), abs(vel_adv_z));
      if delta_t>dt_max then delta_t = dt_max;
      if stepNumber%freq==0 || stepNumber==max_steps then 
        writeln("time ", time, " (step ", stepNumber ,"): dt = ", delta_t);
      if shouldOutput(time, output_interval, delta_t) && (stepNumber>1 || max_steps<10000000000000) {
        prims_to_consv(grid);
        dump_ascii_to_disk(grid, stepNumber, time);
      }
    }
    /* end of computation */
    if max_steps>0 then dump_ascii_to_disk(grid, stepNumber, time);

}
