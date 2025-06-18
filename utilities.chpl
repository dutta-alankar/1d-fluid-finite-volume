module utilities {

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
}