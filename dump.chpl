module dump {
    use globals;
    use structures;
    use utilities;
    use Math;

    var ascii_dump_counter: int(64) = 0;
    var ascii_last_dump_step: int(64) = -1;

    proc dump_ascii_to_disk (grid: borrowed Grid(?), stepNumber: int(64), time: real(64)): void {
        if ascii_last_dump_step==stepNumber then return;
        var computeDomain: domain(3) = {grid.indicesInner.low[0]..grid.indicesInner.high[0], grid.indicesInner.low[1]..grid.indicesInner.high[1], grid.indicesInner.low[2]..grid.indicesInner.high[2]};
        var num_cells = computeDomain.size;
        var positions_x: [0..num_cells-1] real(64);
        var positions_y: [0..num_cells-1] real(64);
        var positions_z: [0..num_cells-1] real(64);
        var state: [0..num_cells-1] real(64);
        writeln("time ", time, " (step ", stepNumber, "): Dumping ./output.", padWithZeros(ascii_dump_counter, 4), ".txt to disk");
        var idx: int = 0;
        for (i,j,k) in computeDomain {
            positions_x[idx] = grid.cells_tot[i,j,k].center[0];
            positions_y[idx] = grid.cells_tot[i,j,k].center[1];
            positions_z[idx] = grid.cells_tot[i,j,k].center[2];
            state[idx] = grid.cells_tot[i,j,k].state_prims_center[grid.cells_tot[i,j,k].states_count.low];
            idx += 1;
        }
        // For now, dumping only x-positions and state. This needs to be generalized for 3D visualization.
        writeArraysToFile("./output."+padWithZeros(ascii_dump_counter, 4)+".txt", positions_x, state);
        ascii_dump_counter += 1;
        ascii_last_dump_step = stepNumber;
    }

    proc shouldOutput (time: real(64), output_dt: real(64), dt: real(64)): bool {
        var time_start: real(64) = 0.0;
        // return time>=(floor((time-time_start)/output_dt)*output_dt) && (time+dt)<(floor((time-time_start)/output_dt)*output_dt);
        return abs(mod(time, output_dt)) < 1.0001*dt;
    }
}