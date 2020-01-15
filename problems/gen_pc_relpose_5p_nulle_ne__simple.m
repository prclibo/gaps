opt = default_options();
opt.M2_path = '/Users/li/workspace/Macaulay2-1.13/bin/M2';
opt.optimize_coefficients = true;

prob_fn = @prob_pc_relpose_5p_nulle_ne__simple;
[solv, opt] = generate_solver(prob_fn, opt);