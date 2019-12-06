opt = default_options();
opt.M2_path = '/Users/li/workspace/Macaulay2-1.13/bin/M2';
% opt.use_sym = false;
% opt.find_sym = false;
% opt.fast_monomial_extraction = false;
% opt.remove_zero_sol = true;
opt.optimize_coefficients = true;

prob_fn = @prob_pc_relpose_4pra_p6d;
prob_fn = @prob_pc_relpose_5p_nulle_ne;
[solv, opt] = generate_solver(prob_fn, opt);

addpath solvers

solv_fun = str2func(['solver_' solv.name]);
stats = benchmark_solver(solv_fun,solv,500);

figure(1)
clf
hist(log10(stats.all_res),50)
title(sprintf('Mean: %.2f, Median: %.2f\n Mode: %.2f, Time: %.2f ms\n',...
        stats.res_mean,stats.res_median,stats.res_mode,...
        1000*median(stats.time_taken)));
xlabel('log10 residual')
ylabel('freq.')
