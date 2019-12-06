function str = cg_compute_coeff(solv, template, opt)

prob = solv.prob;
in_par_str = strjoin(fieldnames(prob.in_subs), ', ');
str = sprintf('function coeffs = compute_coeffs(%s)\n', in_par_str);

par_names = fieldnames(prob.in_subs);
for i = 1:numel(par_names)
    par_name = par_names{i};
    subs_ = prob.in_subs.(par_name);
    for j = 1:numel(subs_)
        line = sprintf('%s = %s(%d);\n', char(subs_(j)), par_name, j);
        str = [str, line];
    end
end

coeff_eqs = solv.coefficients.coeff_eqs;
if opt.optimize_coefficients
    str = [str cg_optimize_coefficients(coeff_eqs)];
else
    expd_str = arrayfun(@char, coeff_eqs, 'UniformOutput', false);
    idx_str = strseq('', 1:numel(coeff_eqs));
    coeff_str = strcat('coeffs(', idx_str(:), ') =\t', expd_str(:), ';\n');
    
    for kk = nvars(coeff_eqs(1)):-1:1
        var_name = sprintf('x%d',kk);
        rep_name = sprintf('data(%d)',kk);
        coeff_str = strrep(coeff_str,var_name,rep_name);
    end
    str = [str coeff_str];
end

end

