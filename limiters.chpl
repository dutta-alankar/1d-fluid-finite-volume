module limiters {
    use Math;

    proc minmod (a: real(64), b: real(64)): real(64) {
        if abs(a)<abs(b) && a*b>0 then
            return a;
        else if abs(b)<abs(a) && a*b>0 then
            return b;
        else
            return 0.0;
    }
    
    proc mc_lim (a: real(64), b: real(64), c: real(64)): real(64) {
        // XXX: sgn might be deprecated later
        var xi: real(64) = b*c;
        if xi > 0 then
            return min(abs(a), abs(b), abs(c))*sgn(a);
        else
            return 0.0;
    }
}