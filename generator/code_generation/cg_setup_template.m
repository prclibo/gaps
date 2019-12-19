function str = cg_setup_template(solv, templ, opt)
prob = solv.prob;
in_par_str = strjoin(fieldnames(prob.in_subs), ', ');
str = sprintf('function [C0,C1] = setup_elimination_template(%s)\n', in_par_str);
str = [str sprintf('[coeffs] = compute_coeffs(%s);\n', in_par_str)];

    function str = cg_setup_block(r, c)
        blk_cc = templ.cind_blk{r, c};
        blk_ind = find(blk_cc);
        blk_ci = blk_cc(blk_ind);
        str = '';
        nm_ind = sprintf('A%d%d_ind', r, c);
        str = [str, cg_format_vector(nm_ind, blk_ind, opt.max_str_len)];
        nm_ci = sprintf('A%d%d_ci', r, c);
        str = [str, cg_format_vector(nm_ci, blk_ci, opt.max_str_len)];
        nm = sprintf('A%d%d', r, c);
        str = [str, sprintf('%s = zeros(%d, %d);\n', nm, size(blk_cc))];
        str = [str, sprintf('%s(%s) = coeffs(%s);\n', nm, nm_ind, nm_ci)];
    end

if opt.find_upper_trianglar
    for r = 1:2
        for c = 1:3
            str = [str, cg_setup_block(r, c)];
        end
    end
    str = [str, sprintf('mplier = A21 / A11;\n')];
    str = [str, sprintf('C0 = A22 - mplier * A12;\n')];
    str = [str, sprintf('C1 = A23 - mplier * A13;\n')];
else
    str = [str cg_format_vector('coeffs0_ind',templ.C0_coeff,opt.max_str_len)];
    str = [str cg_format_vector('coeffs1_ind',templ.C1_coeff,opt.max_str_len)];
    if opt.sparse_template
        str = [str cg_format_vector('ii0',templ.C0_ind2(:,1),opt.max_str_len)];
        str = [str cg_format_vector('jj0',templ.C0_ind2(:,2),opt.max_str_len)];
        str = [str cg_format_vector('ii1',templ.C1_ind2(:,1),opt.max_str_len)];
        str = [str cg_format_vector('jj1',templ.C1_ind2(:,2),opt.max_str_len)];
        str = [str sprintf('C0 = sparse(ii0,jj0,coeffs(coeffs0_ind),%d,%d);\n',templ.C0_sz(1),templ.C0_sz(2))];
        str = [str sprintf('C1 = sparse(ii1,jj1,coeffs(coeffs1_ind),%d,%d);\n',templ.C1_sz(1),templ.C1_sz(2))];
        
    else
        str = [str cg_format_vector('C0_ind',templ.C0_ind,opt.max_str_len)];
        str = [str cg_format_vector('C1_ind',templ.C1_ind,opt.max_str_len)];
        str = [str sprintf('C0 = zeros(%d,%d);\n',templ.C0_sz(1),templ.C0_sz(2))];
        str = [str sprintf('C1 = zeros(%d,%d);\n',templ.C1_sz(1),templ.C1_sz(2))];
        str = [str sprintf('C0(C0_ind) = coeffs(coeffs0_ind);\n')];
        str = [str sprintf('C1(C1_ind) = coeffs(coeffs1_ind);\n')];
        
    end
end

end