opt = default_options();
opt.M2_path = '/Users/li/workspace/Macaulay2-1.13/bin/M2';
opt.optimize_coefficients = true;
% opt.use_sym = false;
opt.force_vars_in_reducibles = true;
prob_fn = @prob_axb_seq_prod;
[solv, opt] = generate_solver(prob_fn, opt);