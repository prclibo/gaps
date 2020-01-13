function str = cg_cpp_setup_templ(solv, templ, opt)
prob = solv.prob;

    function str = cg_cpp_setup_block(r, c)
        blk_cc = templ.cind_blk{r, c};
        blk_ind = find(blk_cc);
        blk_ci = blk_cc(blk_ind) - 1;
        blk_ind = blk_ind - 1;
        str = '';
        nm_ind = sprintf('A%d%d_ind', r, c);
        str = [str, cg_cpp_format_vector(nm_ind, blk_ind, opt.max_str_len)];
        nm_ci = sprintf('A%d%d_ci', r, c);
        str = [str, cg_cpp_format_vector(nm_ci, blk_ci, opt.max_str_len)];
        nm = sprintf('A%d%d', r, c);
        str = [str, sprintf('Eigen::MatrixXd %s = Eigen::MatrixXd::Zero(%d, %d);\n',...
            nm, size(blk_cc))];
        str = [str, sprintf('for (int i = 0; i < %d; ++i) %s(%s[i]) = coeffs(%s[i]);\n',...
            numel(blk_ind), nm, nm_ind, nm_ci)];
    end
str = '';
if opt.find_upper_trianglar
    for r = 1:2
        for c = 1:3
            str = [str, cg_cpp_setup_block(r, c)];
        end
    end
    str = [str, sprintf('Eigen::MatrixXd mplier = A11.transpose().fullPivLu().solve(A21.transpose()).transpose();\n')];
    str = [str, sprintf('Eigen::MatrixXd C0 = A22 - mplier * A12;\n')];
    str = [str, sprintf('Eigen::MatrixXd C1 = A23 - mplier * A13;\n')];
    str = [str, sprintf('MatrixXd C12 = C0.fullPivLu().solve(C1);\n')];
    
else
    if opt.sparse_template
        error('libo: Not fixed yet');
        % Make sure we sort by columns since the sparse matrix is stored in
        % column-major order.
        [outerInd,ind] = sort(templ.C0_ind2(:,2));
        innerInd = templ.C0_ind2(ind,1);
        outerStarts = cumsum([1 histcounts(outerInd,(1:templ.C0_sz(2)+1)-0.5)]);
        
        str = [str cg_cpp_format_vector('C0_outer_indices',outerStarts-1,opt.max_str_len) sprintf('\n')];
        str = [str cg_cpp_format_vector('C0_inner_indices',innerInd-1,opt.max_str_len) sprintf('\n')];
        str = [str sprintf('VectorXd C0_values(%d);\n',length(templ.C0_ind))];
        str = [str sprintf('for (int i = 0; i < %d; i++) {\n',length(templ.C0_ind))];
        str = [str opt.cg_indentation sprintf('C0_values[i] = coeffs(coeffs0_ind[i]);\n}\n')];
        str = [str sprintf('const SparseMatrix<double> C0 = Map<const SparseMatrix<double>>(%d, %d, %d, C0_outer_indices, C0_inner_indices, C0_values.data());\n\n',...
            templ.C0_sz(1),templ.C0_sz(2),length(templ.C0_ind))];
        
        [outerInd,ind] = sort(templ.C1_ind2(:,2));
        innerInd = templ.C1_ind2(ind,1);
        outerStarts = cumsum([1 histcounts(outerInd,(1:templ.C1_sz(2)+1)-0.5)]);
        
        str = [str cg_cpp_format_vector('C1_outer_indices',outerStarts-1,opt.max_str_len) sprintf('\n')];
        str = [str cg_cpp_format_vector('C1_inner_indices',innerInd-1,opt.max_str_len) sprintf('\n')];
        str = [str sprintf('VectorXd C1_values(%d);\n',length(templ.C1_ind))];
        str = [str sprintf('for (int i = 0; i < %d; i++) {\n',length(templ.C1_ind))];
        str = [str opt.cg_indentation sprintf('C1_values[i] = coeffs(coeffs0_ind[i]);\n}\n')];
        str = [str sprintf('const SparseMatrix<double> C1 = Map<const SparseMatrix<double>>(%d, %d, %d, C1_outer_indices, C1_inner_indices, C1_values.data());\n\n',...
            templ.C1_sz(1),templ.C1_sz(2),length(templ.C1_ind))];
        % Construct C from triplets instead of mapping to data.
        %     str = [];
        %     str = [str 'Triplet<double> C_triplets[] = {' sprintf('{%d,%d,%d},',[innerInd'-1;outerInd'-1;-ones(1,length(innerInd))]) '};' newline];
        %     str = [str sprintf('SparseMatrix<double> C(%d,%d);\n',templ.C_sz(1),templ.C_sz(2))];
        %     str = [str sprintf('C.setFromTriplets(&C_triplets[0], &C_triplets[%d]);\n\n',length(innerInd))];
    else
        str = [str cg_cpp_format_vector('C0_ind',templ.C0_ind-1,opt.max_str_len) sprintf('\n')];
        str = [str cg_cpp_format_vector('C1_ind',templ.C1_ind-1,opt.max_str_len) sprintf('\n')];
        str = [str sprintf('MatrixXd C0 = MatrixXd::Zero(%d,%d);\n',templ.C0_sz(1),templ.C0_sz(2))];
        str = [str sprintf('MatrixXd C1 = MatrixXd::Zero(%d,%d);\n',templ.C1_sz(1),templ.C1_sz(2))];
        str = [str sprintf('for (int i = 0; i < %d; i++) {\n',length(templ.C0_ind))];
        str = [str opt.cg_indentation sprintf('C0(C0_ind[i]) = coeffs(coeffs0_ind[i]);\n}\n\n')];
        str = [str sprintf('for (int i = 0; i < %d; i++) {\n',length(templ.C1_ind))];
        str = [str opt.cg_indentation sprintf('C1(C1_ind[i]) = coeffs(coeffs1_ind[i]);\n}\n\n')];
        str = [str, sprintf('MatrixXd C12 = C0.fullPivLu().solve(C1);\n')];

    end
    
end

end