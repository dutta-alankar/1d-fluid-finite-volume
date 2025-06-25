module dump {
    use structures;
    use utilities;
    use Math;

    var ascii_dump_counter: int(64) = 0;
    var ascii_last_dump_step: int(64) = -1;

    proc dump_ascii_to_disk (grid: borrowed Grid(?), stepNumber: int(64), time: real(64)): void {
        if ascii_last_dump_step==stepNumber then return;
        var computeDomain: domain(1) = {grid.indicesInner.low..grid.indicesInner.high};
        var positions: [computeDomain] real(64);
        var state: [computeDomain] real(64);
        writeln("time ", time, " (step ", stepNumber, "): Dumping ./output.", padWithZeros(ascii_dump_counter, 4), ".txt to disk");
        forall i in computeDomain {
            positions[i] = grid.cells_tot[i].center;
            state[i] = grid.cells_tot[i].state_prims_center[grid.cells_tot[i].states_count.low];
        }
        writeArraysToFile("./output."+padWithZeros(ascii_dump_counter, 4)+".txt", positions, state);
        ascii_dump_counter += 1;
        ascii_last_dump_step = stepNumber;
    }

    proc shouldOutput (time: real(64), output_dt: real(64), dt: real(64)): bool {
        return abs(mod(time, output_dt)) < 1.1*dt;
    }
}