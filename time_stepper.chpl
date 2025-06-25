module time_stepper {
    use structures;

    proc timeStep (consv_state: [state_domain] real(64), rhs_state: [state_domain] real(64), dt: real(64)): [state_domain] real(64) {
        /* explicit Euler first order time integration */
        var consv_state_final: [state_domain] real(64);
        for state_var in state_domain {
            consv_state_final[state_var] = consv_state[state_var] + rhs_state[state_var]*dt;
        }
        return consv_state_final;
    }
}