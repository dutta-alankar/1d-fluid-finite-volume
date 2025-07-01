module globals {
    config const npoints_x: int(64)  = 5;
    config const npoints_y: int(64)  = 1;
    config const npoints_z: int(64)  = 1;
    config const xmin_x:    real(64) = 0.0;
    config const xmin_y:    real(64) = 0.0;
    config const xmin_z:    real(64) = 0.0;
    config const xmax_x:    real(64) = 5.0;
    config const xmax_y:    real(64) = 1.0;
    config const xmax_z:    real(64) = 1.0;
    config const reconstruction_type: string = "constant";
    config const time_integration: string = "euler";
    config const cfl: real(64) = 0.3;
    config const dt_ini: real(64) = 1.0e-04;
    config const t_start: real(64) = 0.0;
    config const t_stop: real(64) = 3.0;
    config const vel_adv_x: real(64) = 1.0;
    config const vel_adv_y: real(64) = 0.0;
    config const vel_adv_z: real(64) = 0.0;
    config const dt_max: real(64) = 1.0e-02;
    config const max_steps: int(64) = 10000000000000;
    config const freq: int(64) = 20;
    config const output_interval: real(64) = 0.1;
    config const slope_limiter: string = "mc";
}