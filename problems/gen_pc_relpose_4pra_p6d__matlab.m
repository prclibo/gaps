opt = default_options();
opt.M2_path = '/Users/li/workspace/Macaulay2-1.13/bin/M2';
opt.optimize_coefficients = true;
opt.remove_extra_columns = false;
opt.find_upper_trianglar = true;

opt.find_upper_trianglar = false;

prob_fn = @prob_pc_relpose_4pra_p6d;
[solv, opt] = generate_solver(prob_fn, opt);
% return;
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
