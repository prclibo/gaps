function str = cg_cpp_setup_template(solv, template, opt)


str = [cg_cpp_format_vector('coeffs0_ind',template.C0_coeff-1,opt.max_str_len)];
str = [str sprintf('\n')];
str = [str cg_cpp_format_vector('coeffs1_ind',template.C1_coeff-1,opt.max_str_len)];
str = [str sprintf('\n')];

if opt.sparse_template
    % Make sure we sort by columns since the sparse matrix is stored in
    % column-major order.
    [outerInd,ind] = sort(template.C0_ind2(:,2));
    innerInd = template.C0_ind2(ind,1);
    outerStarts = cumsum([1 histcounts(outerInd,(1:template.C0_sz(2)+1)-0.5)]);

    str = [str cg_cpp_format_vector('C0_outer_indices',outerStarts-1,opt.max_str_len) sprintf('\n')];
    str = [str cg_cpp_format_vector('C0_inner_indices',innerInd-1,opt.max_str_len) sprintf('\n')];
    str = [str sprintf('VectorXd C0_values(%d);\n',length(template.C0_ind))];
    str = [str sprintf('for (int i = 0; i < %d; i++) {\n',length(template.C0_ind))];
    str = [str opt.cg_indentation sprintf('C0_values[i] = coeffs(coeffs0_ind[i]);\n}\n')];
    str = [str sprintf('const SparseMatrix<double> C0 = Map<const SparseMatrix<double>>(%d, %d, %d, C0_outer_indices, C0_inner_indices, C0_values.data());\n\n',...
        template.C0_sz(1),template.C0_sz(2),length(template.C0_ind))];
    
    [outerInd,ind] = sort(template.C1_ind2(:,2));
    innerInd = template.C1_ind2(ind,1);
    outerStarts = cumsum([1 histcounts(outerInd,(1:template.C1_sz(2)+1)-0.5)]);

    str = [str cg_cpp_format_vector('C1_outer_indices',outerStarts-1,opt.max_str_len) sprintf('\n')];
    str = [str cg_cpp_format_vector('C1_inner_indices',innerInd-1,opt.max_str_len) sprintf('\n')];
    str = [str sprintf('VectorXd C1_values(%d);\n',length(template.C1_ind))];
    str = [str sprintf('for (int i = 0; i < %d; i++) {\n',length(template.C1_ind))];
    str = [str opt.cg_indentation sprintf('C1_values[i] = coeffs(coeffs0_ind[i]);\n}\n')];
    str = [str sprintf('const SparseMatrix<double> C1 = Map<const SparseMatrix<double>>(%d, %d, %d, C1_outer_indices, C1_inner_indices, C1_values.data());\n\n',...
        template.C1_sz(1),template.C1_sz(2),length(template.C1_ind))];
    % Construct C from triplets instead of mapping to data.
%     str = [];
%     str = [str 'Triplet<double> C_triplets[] = {' sprintf('{%d,%d,%d},',[innerInd'-1;outerInd'-1;-ones(1,length(innerInd))]) '};' newline];
%     str = [str sprintf('SparseMatrix<double> C(%d,%d);\n',template.C_sz(1),template.C_sz(2))];
%     str = [str sprintf('C.setFromTriplets(&C_triplets[0], &C_triplets[%d]);\n\n',length(innerInd))];
else
    str = [str cg_cpp_format_vector('C0_ind',template.C0_ind-1,opt.max_str_len) sprintf('\n')];
    str = [str cg_cpp_format_vector('C1_ind',template.C1_ind-1,opt.max_str_len) sprintf('\n')];
    str = [str sprintf('MatrixXd C0 = MatrixXd::Zero(%d,%d);\n',template.C0_sz(1),template.C0_sz(2))];
    str = [str sprintf('MatrixXd C1 = MatrixXd::Zero(%d,%d);\n',template.C1_sz(1),template.C1_sz(2))];
    str = [str sprintf('for (int i = 0; i < %d; i++) {\n',length(template.C0_ind))];
    str = [str opt.cg_indentation sprintf('C0(C0_ind[i]) = coeffs(coeffs0_ind[i]);\n}\n\n')];
    str = [str sprintf('for (int i = 0; i < %d; i++) {\n',length(template.C1_ind))];
    str = [str opt.cg_indentation sprintf('C1(C1_ind[i]) = coeffs(coeffs1_ind[i]);\n}\n\n')];    
end

