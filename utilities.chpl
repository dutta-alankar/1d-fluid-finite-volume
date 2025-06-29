module utilities {
    use IO;
    use Math;

    proc linspace(start: real(64), stop: real(64), num: int(64), D: domain(?)): [D] real(64) {
        assert(D.rank == 1, "Domain must be one-dimensional 'rank'="+D.rank:string+"!=1");
        assert(D.size == num, "Domain size "+D.size:string+" must match 'num'="+num:string);
        var result: [D] real(64);
        if (num == 1) then {
            result[D.low] = start;
        } else {
            const step = (stop - start) / (num - 1): real(64);
            forall i in D {
                result[i] = start + (i - D.low) * step;
            }
        }
        return result;
    }

    // Parallel write using zippered iteration
    proc writeArraysToFile(filename: string, A: [] real(64), B: [] real(64)) {
        assert(A.domain==B.domain && A.domain.rank==1, "Domains of A and B must match and be one-dimensional. A:"+A.domain:string+", B:"+B.domain:string);
        var f = try! open(filename, ioMode.cw);
        var writer = try! f.writer(locking=true);
        coforall loc in Locales {
            on loc {
                forall (a, b, i) in zip(A, B, A.domain) with (ref writer) {
                    try! writer.writeln(a, "\t", b);
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