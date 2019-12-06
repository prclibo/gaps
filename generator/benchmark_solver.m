function [ result ] = benchmark_solver( solv_fun, solv, iters, discard_zero_sol )

if nargin < 3
    iters = 10;
end

if nargin < 4
    discard_zero_sol = 1;
end

result = [];

result.all_res = [];
result.failures = 0;
result.time_taken = [];

% Hack to figure out how large the data is
% TODO make better choices in design (and life)
prob = solv.prob;

fprintf('Solving sample instances for           ');
for iter = 1:iters
    fprintf('\b\b\b\b\b\b\b\b%8d', iter);

    rng(iter);
    [in_rl, out_rl] = prob.rand_par_rl();

    out_cell = cell(size(fieldnames(out_rl)));
    in_cell = struct2cell(in_rl);
    
    try
        tic;
        [out_cell{:}] = solv_fun(in_cell{:});
        tt = toc;
    catch
        result.failures = result.failures + 1;
        continue;
    end

    result.time_taken(end+1) = tt;    

%     if discard_zero_sol
%         sols = sols(:,max(abs(sols))>1e-10);
%     end
    
%     % We measure maximum equation residual
%     res = [];
%     for k = 1:size(sols,2)
%         res = [res max(abs(evaluate(eqs,sols(:,k))))];
%     end

    % Measure minimal distance to groundtruth solution
    min_d = inf;
    for k = 1:numel(out_cell{1})
        sol = cellfun(@(x) x{k}, out_cell, 'UniformOutput', false);
        sol = cat(1, sol{:});
        gt = struct2cell(out_rl);
        gt = cat(1, gt{:});
        d = norm(sol - gt);
        min_d = min(min_d, d);
    end
    
    result.all_res = [result.all_res min_d];
end
fprintf('\n');

result.res_mean = mean(log10(result.all_res));
result.res_median = median(log10(result.all_res));
[hh,bb]=hist(log10(result.all_res),20);
[~,idde]=max(hh);
result.res_mode = bb(idde);

end

