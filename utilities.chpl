module utilities {
    use IO;
    use Math;

    // Parallel write using zippered iteration
    proc writeArraysToFile(filename: string, A: [] real(64), B: [] real(64)) {
        assert(A.domain==B.domain && A.domain.rank==1, "Domains of A and B must match and be one-dimensional. A:"+A.domain:string+", B:"+B.domain:string);
        var f = try! open(filename, ioMode.cw);
        var writer = try! f.writer(locking=true);
        coforall loc in Locales {
            on loc {
                forall (a, b, i) in zip(A, B, A.domain) with (ref writer) {
                    try! writer.writeln(a, "	", b);
                }
            }
        }
        sync try! writer.close();
    }

    proc padWithZeros(n: int(64), width: int(64)): string {
        var s: string = n:string;
        var zerosNeeded: int(64) = width - s.numBytes;
        if zerosNeeded > 0 then
            return "0" * zerosNeeded + s;
        else 
            return s;
    }

}