clear
opt = default_options();
opt.M2_path = '/Users/li/workspace/Macaulay2-1.13/bin/M2';
opt.optimize_coefficients = true;
opt.use_sym = false;
opt.force_vars_in_reducibles = true;
opt.integer_expansions = 6;
opt.remove_extra_columns = false;
% opt.find_upper_trianglar = true;

config.rand_seed = 2;

opt.solver_name = sprintf('xinvax_eq_b');

prob_fn = @() prob_xinvax_eq_b(config);
[solv, opt] = generate_solver(prob_fn, opt);

% return;
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
