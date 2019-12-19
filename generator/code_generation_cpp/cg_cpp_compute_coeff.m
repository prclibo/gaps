function str = cg_cpp_compute_coeff(solv, template, opt)

prob = solv.prob;

prob = solv.prob;
in_par_str = strjoin(fieldnames(prob.in_subs), ', ');
str = '';

par_names = fieldnames(prob.in_subs);
for i = 1:numel(par_names)
    par_name = par_names{i};
    subs_ = prob.in_subs.(par_name);
    if numel(subs_) > 1
        str = [str, sprintf('double const* %s_data = %s.data();\n', par_name, par_name)];
        for j = 1:numel(subs_)
            line = sprintf('double const %s = %s_data[%d];\n', char(subs_(j)), par_name, j - 1);
            str = [str, line];
        end
    end
end

coeff_eqs = solv.coefficients.coeff_eqs;
if opt.optimize_coefficients
    str = [str sprintf('%sVectorXd coeffs(%d);\n',opt.cg_indentation,length(coeff_eqs))];
    str = [str cg_cpp_optimize_coefficients(coeff_eqs)];
else
    % Accessing the underlying data directly will significantly
    % improve the compilation time, rather than accessing the
    % elements of the VectorXd data.
    str = [str sprintf('%sconst double* d = data.data();\n',opt.cg_indentation)];
    str = [str sprintf('%sVectorXd coeffs(%d);\n',opt.cg_indentation,length(coeff_eqs))];
    coeff_str = [];
    for k = 1:length(coeff_eqs)
        c = char(coeff_eqs(k),0);
        coeff_str = [coeff_str sprintf('%scoeffs[%d] = %s;\n',opt.cg_indentation,k-1,c)];
    end
    
    for kk = nvars(coeff_eqs(1)):-1:1
        var_name = sprintf('x%d',kk);
        rep_name = sprintf('d[%d]',kk-1);
        coeff_str = strrep(coeff_str,var_name,rep_name);
    end
    
    % The compiler will optimize away the std::pow calls for low
    % exponents, eg. std::pow(x,2) -> x*x.
    coeff_str = regexprep(coeff_str,'(d\[\d+\])\^(\d+)','std::pow($1,$2)');
    
    str = [str coeff_str sprintf('\n')];
end
end

