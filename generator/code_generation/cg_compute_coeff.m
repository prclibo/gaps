function str = cg_compute_coeff(solv, template, opt)

prob = solv.prob;
in_par_str = strjoin(fieldnames(prob.in_subs), ', ');
str = sprintf('function coeffs = compute_coeffs(%s)\n', in_par_str);

% Unpack input parameters
par_names = fieldnames(prob.in_subs);
for i = 1:numel(par_names)
    par_name = par_names{i};
    subs_ = prob.in_subs.(par_name);
    for j = 1:numel(subs_)
        line = sprintf('%s = %s(%d);\n', char(subs_(j)), par_name, j);
        str = [str, line];
    end
end

par_names = fieldnames(prob.abbr_subs);
for i = 1:numel(par_names)
    par_name = par_names{i};
    subs_ = prob.abbr_subs.(par_name);
    line = sprintf('%s = %s;\n', par_name, char(subs_));
    str = [str, line];
end

coeff_eqs = solv.coefficients.coeff_eqs;
if opt.optimize_coefficients
    str = [str cg_optimize_coefficients(coeff_eqs)];
else
    expd_str = arrayfun(@char, coeff_eqs, 'UniformOutput', false);
    idx_str = strseq('', 1:numel(coeff_eqs));
    coeff_str = strcat('coeffs(', idx_str(:), ') =\t', expd_str(:), ';\n');
    coeff_str = sprintf(horzcat(coeff_str{:}));
    str = [str coeff_str];
end

end

