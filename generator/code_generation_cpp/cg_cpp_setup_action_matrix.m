function str = cg_cpp_setup_action_matrix(solv,template,opt)
str = [];

if opt.generalized_eigenvalue_solver
    % TODO: Implement when no longer experimental feature.
    error('Generalized eigenvalue solver has not been implemented for C++ output.');
end

if template.C_sz(1) == template.C_sz(2)-length(template.basis)
    % Square template
    if opt.sparse_template
       
        if strcmp(opt.linear_solver,'backslash') || strcmp(opt.linear_solver,'lu')
            str = [str 'SparseLU<SparseMatrix<double>,COLAMDOrdering<SparseMatrix<double>::StorageIndex>> solver;' sprintf('\n')];
        elseif strcmp(opt.linear_solver,'qr')
            str = [str 'SparseQR<SparseMatrix<double>,COLAMDOrdering<SparseMatrix<double>::StorageIndex>> solver;' sprintf('\n')];
        elseif strcmp(opt.linear_solver,'cg')
            str = [str 'LeastSquaresConjugateGradient<SparseMatrix<double>> solver;' ];
        else
            error([opt.linear_solver ' is an unknown or unimplemented solver for sparse C++ output.']);
        end
        str = [str 'solver.compute(C0);' sprintf('\n')];
        str = [str 'MatrixXd C12 = solver.solve(C1);' sprintf('\n')];
    else        
        if strcmp(opt.linear_solver,'backslash') || strcmp(opt.linear_solver,'lu')
            %str = [str sprintf('MatrixXd C12 = C0.fullPivLu().solve(C1);\n')];
            str = [str sprintf('MatrixXd C12 = C0.partialPivLu().solve(C1);\n')];
        elseif strcmp(opt.linear_solver,'qr')
            str = [str sprintf('MatrixXd C12 = C0.fullPivHouseholderQr().solve(C1);\n')];
        else
            error([opt.linear_solver ' is an unknown or unimplemented solver for dense C++ output.']);
        end
    end
    str = [str sprintf('MatrixXd RR(%d, %d);\n',length(template.reducible)+length(template.basis),length(template.basis))];
    str = [str sprintf('RR << -C12.bottomRows(%d), MatrixXd::Identity(%d, %d);\n\n',length(template.reducible),length(template.basis),length(template.basis))];
else
    b = zeros(template.C_sz(2)-length(template.basis),length(template.reducible));
    b(end-length(template.reducible)+1:end,:) = -eye(length(template.reducible));
    b_ind = find(b);
    b_sz = size(b);

    % Sparse b does not work for the SparseQR solver which is why we aways
    % use a dense representation. However, a sparse b might be beneficial
    % for other solvers.
    str = [str cg_cpp_format_vector('b_ind',b_ind-1,opt.max_str_len) sprintf('\n')];
    str = [str sprintf('MatrixXd b = MatrixXd::Zero(%d,%d);\n',b_sz(1),b_sz(2))];
    str = [str sprintf('for (int i = 0; i < %d; i++) {\n',length(b_ind))];
    str = [str opt.cg_indentation sprintf('b(b_ind[i]) = -1;\n}\n\n')];

    if opt.sparse_template
        % Construct b as a sparse matrix.
%         [innerInd,outerInd] = ind2sub(b_sz,b_ind);
%         outerStarts = cumsum([1 histcounts(outerInd,(1:b_sz(2)+1)-0.5)]);
%         str = [str cg_cpp_format_vector('b_outer_indices',outerStarts-1,opt.max_str_len)];
%         str = [str cg_cpp_format_vector('b_inner_indices',innerInd-1,opt.max_str_len)];
%         str = [str 'static const double b_values[] = {' sprintf('%d,',-ones(1,length(b_ind)-1)) '-1};' sprintf('\n')];
%         str = [str sprintf('const SparseMatrix<double> b = Map<const SparseMatrix<double>>(%d, %d, %d, b_outer_indices, b_inner_indices, b_values);\n\n',...
%             b_sz(1),b_sz(2),length(b_ind))];

        if strcmp(opt.linear_solver,'backslash') || strcmp(opt.linear_solver,'qr')
            str = [str 'SparseQR<SparseMatrix<double>,COLAMDOrdering<SparseMatrix<double>::StorageIndex>> solver;' sprintf('\n')];
        elseif strcmp(opt.linear_solver,'cg')
            str = [str 'LeastSquaresConjugateGradient<SparseMatrix<double>> solver;' sprintf('\n')];
        else
            error([opt.linear_solver ' is an unknown or unimplemented solver for sparse C++ output.']);
        end
        str = [str 'solver.compute(C0.adjoint());' sprintf('\n')];
        str = [str 'MatrixXd alpha = solver.solve(b);' sprintf('\n')];
    else
        if strcmp(opt.linear_solver,'backslash') || strcmp(opt.linear_solver,'lu')
            % TODO: For this to be equivalent to Matlabs "backslash"
            % operator there is a lot of logic that is needed here. For
            % example, C0 might be square, symmetric, SPD, triangular etc.
            str = [str sprintf('MatrixXd alpha = C0.adjoint().fullPivLu().solve(b);\n')];
        elseif strcmp(opt.linear_solver,'qr')
            str = [str sprintf('MatrixXd alpha = C0.adjoint().fullPivHouseholderQr().solve(b);\n')];
        elseif strcmp(opt.linear_solver,'svd')
            str = [str sprintf('MatrixXd alpha = C0.adjoint().jacobiSvd(ComputeThinU | ComputeThinV).solve(b);\n')];
        else
            error([opt.linear_solver ' is an unknown or unimplemented solver for dense C++ output.']);
        end
    end
    str = [str sprintf('MatrixXd RR(%d, %d);\n',length(b_ind)+length(template.basis),length(template.basis))];
    str = [str sprintf('RR << alpha.adjoint()*C1, MatrixXd::Identity(%d, %d);\n\n',length(template.basis),length(template.basis))];
end

str = [str cg_cpp_format_vector('AM_ind',template.AM_ind-1,opt.max_str_len)];    
str = [str sprintf('MatrixXd AM(%d, %d);\n',length(template.AM_ind),length(template.basis))];
str = [str sprintf('for (int i = 0; i < %d; i++) {\n',length(template.AM_ind))];
str = [str opt.cg_indentation sprintf('AM.row(i) = RR.row(AM_ind[i]);\n}\n\n')];
end


