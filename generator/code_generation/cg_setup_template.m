function str = cg_setup_template(solv, template, opt)
prob = solv.prob;
in_par_str = strjoin(fieldnames(prob.in_subs), ', ');
str = sprintf('function [C0,C1] = setup_elimination_template(%s)\n', in_par_str);
str = [str sprintf('[coeffs] = compute_coeffs(%s);\n', in_par_str)];

str = [str cg_format_vector('coeffs0_ind',template.C0_coeff,opt.max_str_len)];
str = [str cg_format_vector('coeffs1_ind',template.C1_coeff,opt.max_str_len)];


if opt.sparse_template        
    str = [str cg_format_vector('ii0',template.C0_ind2(:,1),opt.max_str_len)];
    str = [str cg_format_vector('jj0',template.C0_ind2(:,2),opt.max_str_len)];
    str = [str cg_format_vector('ii1',template.C1_ind2(:,1),opt.max_str_len)];
    str = [str cg_format_vector('jj1',template.C1_ind2(:,2),opt.max_str_len)];
    str = [str sprintf('C0 = sparse(ii0,jj0,coeffs(coeffs0_ind),%d,%d);\n',template.C0_sz(1),template.C0_sz(2))];
    str = [str sprintf('C1 = sparse(ii1,jj1,coeffs(coeffs1_ind),%d,%d);\n',template.C1_sz(1),template.C1_sz(2))];
    
else
    str = [str cg_format_vector('C0_ind',template.C0_ind,opt.max_str_len)];
    str = [str cg_format_vector('C1_ind',template.C1_ind,opt.max_str_len)];
    str = [str sprintf('C0 = zeros(%d,%d);\n',template.C0_sz(1),template.C0_sz(2))];
    str = [str sprintf('C1 = zeros(%d,%d);\n',template.C1_sz(1),template.C1_sz(2))];
    str = [str sprintf('C0(C0_ind) = coeffs(coeffs0_ind);\n')];
    str = [str sprintf('C1(C1_ind) = coeffs(coeffs1_ind);\n')];
end
