function template = generate_code_cpp(solv_name,solv,template,opt)

if ~isfield(template,'C_sz') || isempty(template.C_sz)
    fprintf('Template indicies are missing. Generating...\n');
    template = finalize_template(solv,template);
end

addpath([fileparts(mfilename('fullpath')) filesep 'code_generation_cpp'])

template_dir = [fileparts(mfilename('fullpath')) filesep 'code_generation_cpp' filesep 'templates'];

template_file = [template_dir filesep 'template_solver.cpp'];

code_template_str = fileread(template_file);

replacements = {};

reduced_eigenvector_solver = 0;

% Deal with conditionals first
if strcmp(opt.eigen_solver,'default')
    code_template_str = regexprep(code_template_str,'#if \$\(use_standard_eigensolver\)\s*\n(.*?)\n#endif','$1');
    code_template_str = regexprep(code_template_str,'#if \$\(use_eigsonly_eigensolver\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if \$\(use_sturm_eigensolver\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_sturm_dani_eigensolver\)\s*\n(.*?)\n#endif','');
end

if strcmp(opt.eigen_solver,'eigs_only')
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_standard_eigensolver\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_eigsonly_eigensolver\)\s*\n(.*?)\n#endif','$1');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_sturm_eigensolver\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_sturm_dani_eigensolver\)\s*\n(.*?)\n#endif','');
    reduced_eigenvector_solver = 1;
end

if strcmp(opt.eigen_solver,'sturm')
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_standard_eigensolver\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_eigsonly_eigensolver\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_sturm_eigensolver\)\s*\n(.*?)\n#endif','$1');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_sturm_dani_eigensolver\)\s*\n(.*?)\n#endif','');
    reduced_eigenvector_solver = 1;
    
    replacements{end+1} = {'charpoly_method',                 opt.charpoly_method};
end

if strcmp(opt.eigen_solver,'sturm_dani')
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_standard_eigensolver\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_eigsonly_eigensolver\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_sturm_eigensolver\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_sturm_dani_eigensolver\)\s*\n(.*?)\n#endif','$1');
end


if opt.sparse_template
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_dense_template_stack_alloc\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_dense_template_heap_alloc\)\s*\n(.*?)\n#endif','');
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_sparse_template\)\s*\n(.*?)\n#endif','$1');
else
    if max(template.C_sz) > opt.cg_max_stack_alloc_size
        code_template_str = regexprep(code_template_str,'#if\s*\$\(use_dense_template_stack_alloc\)\s*\n(.*?)\n#endif','');
        code_template_str = regexprep(code_template_str,'#if\s*\$\(use_dense_template_heap_alloc\)\s*\n(.*?)\n#endif','$1');
        code_template_str = regexprep(code_template_str,'#if\s*\$\(use_sparse_template\)\s*\n(.*?)\n#endif','');    
    else
        code_template_str = regexprep(code_template_str,'#if\s*\$\(use_dense_template_stack_alloc\)\s*\n(.*?)\n#endif','$1');
        code_template_str = regexprep(code_template_str,'#if\s*\$\(use_dense_template_heap_alloc\)\s*\n(.*?)\n#endif','');
        code_template_str = regexprep(code_template_str,'#if\s*\$\(use_sparse_template\)\s*\n(.*?)\n#endif','');    
    end
end

if reduced_eigenvector_solver
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_reduced_eigenvector_solver\)\s*\n(.*?)\n#endif','$1');

    red_eig_solver = setup_reduced_eigenvector_solver(template);
    
    replacements{end+1} = {'ind_non_trivial',           commavec(red_eig_solver.ind_non_trivial-1)};
    replacements{end+1} = {'length_ind_non_trivial',    length(red_eig_solver.ind_non_trivial)};
    replacements{end+1} = {'AA_sz',                     commavec(red_eig_solver.AA_sz)};
    replacements{end+1} = {'ind_var',                   commavec(red_eig_solver.ind_var-1)};
    replacements{end+1} = {'ind_unit',                  red_eig_solver.ind_unit-1};
    replacements{end+1} = {'max_power',                 red_eig_solver.max_power};

    [code_setup_AA,code_extract_sol] = cg_cpp_reduced_eigenvector_solver(red_eig_solver,opt);
    
    replacements{end+1} = {'code_setup_reduced_eigenvalue_eq',         code_setup_AA};
    replacements{end+1} = {'code_extract_solutions',                   code_extract_sol};

elseif strcmp(opt.eigen_solver,'sturm_dani')
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_reduced_eigenvector_solver\)\s*\n(.*?)\n#endif','');

    % Compute eigenvectors using danilevsky transform
%     replacements{end+1} = {'code_extract_solutions',         cg_cpp_extract_solutions_sturm_dani(solv,template,opt)};
    replacements{end+1} = {'code_extract_solutions',         cg_cpp_extract_solutions(solv,template,opt)};
else
    % Compute eigenvectors using standard approach
    code_template_str = regexprep(code_template_str,'#if\s*\$\(use_reduced_eigenvector_solver\)\s*\n(.*?)\n#endif','');
    replacements{end+1} = {'code_normalize_eigenvectors',    cg_cpp_normalize_eigenvectors(solv,template,opt)};
    replacements{end+1} = {'code_extract_solutions',         cg_cpp_extract_solutions(solv,template,opt)};
end
replacements{end+1} = {'code_pack_outputs',                        cg_cpp_pack_outputs(solv, template, opt)};

% Misc things
replacements{end+1} = {'solv_name',                 solv_name};
replacements{end+1} = {'debug_comments',            print_debug_comments(solv,template,opt)};
replacements{end+1} = {'code_compute_coefficients', cg_cpp_compute_coeff(solv,template,opt)};
replacements{end+1} = {'code_setup_template',       cg_cpp_setup_template(solv,template,opt)};
replacements{end+1} = {'num_vars',                  numel(solv.prob.unk_vars)};
replacements{end+1} = {'num_basis',                 length(template.basis)};
replacements{end+1} = {'num_reducible',             length(template.reducible)};
replacements{end+1} = {'num_available',             length(template.basis)+length(template.reducible)};
replacements{end+1} = {'input_size',                numel(fieldnames(solv.prob.in_subs))};
replacements{end+1} = {'output_size',               numel(fieldnames(solv.prob.out_subs))};


replacements{end+1} = {'function_param_declaration', cg_cpp_declare_function_parameters(solv, true) };
replacements{end+1} = {'function_param',            cg_cpp_declare_function_parameters(solv, false) };
replacements{end+1} = {'code_mex_function_body',    cg_cpp_mex_function_body(solv, opt)};


% Action matrix stuff
replacements{end+1} = {'AM_ind',                    commavec(template.AM_ind-1)};


% Coefficient indices
replacements{end+1} = {'coeffs_ind',                commavec(template.C_coeff-1)};
replacements{end+1} = {'length_coeffs_ind',         length(template.C_coeff)};
replacements{end+1} = {'coeffs0_ind',               commavec(template.C0_coeff-1)};
replacements{end+1} = {'length_coeffs0_ind',        length(template.C0_coeff)};
replacements{end+1} = {'coeffs1_ind',               commavec(template.C1_coeff-1)};
replacements{end+1} = {'length_coeffs1_ind',        length(template.C1_coeff)};


% Elim. template indices
replacements{end+1} = {'length_C_ind',              length(template.C_ind)};
replacements{end+1} = {'C_sz',                      commavec(template.C_sz)};
replacements{end+1} = {'length_C0_ind',             length(template.C0_ind)};
replacements{end+1} = {'C0_sz',                     commavec(template.C0_sz)};
replacements{end+1} = {'length_C1_ind',             length(template.C1_ind)};
replacements{end+1} = {'C1_sz',                     commavec(template.C1_sz)};
if ~opt.sparse_template
    replacements{end+1} = {'C_ind',                     commavec(template.C_ind-1)};   
    replacements{end+1} = {'C0_ind',                    commavec(template.C0_ind-1)};  
    replacements{end+1} = {'C1_ind',                    commavec(template.C1_ind-1)};    
else
    [outerInd,ind] = sort(template.C_ind2(:,2));innerInd = template.C_ind2(ind,1);
    outerStarts = cumsum([1 histcounts(outerInd,(1:template.C_sz(2)+1)-0.5)]);
    replacements{end+1} = {'C_outer_ind',             commavec(outerStarts-1)};
    replacements{end+1} = {'C_inner_ind',             commavec(innerInd-1)};    
    
    [outerInd,ind] = sort(template.C0_ind2(:,2));innerInd = template.C0_ind2(ind,1);
    outerStarts = cumsum([1 histcounts(outerInd,(1:template.C0_sz(2)+1)-0.5)]);
    replacements{end+1} = {'C0_outer_ind',             commavec(outerStarts-1)};
    replacements{end+1} = {'C0_inner_ind',             commavec(innerInd-1)};
    
    [outerInd,ind] = sort(template.C1_ind2(:,2));innerInd = template.C1_ind2(ind,1);
    outerStarts = cumsum([1 histcounts(outerInd,(1:template.C1_sz(2)+1)-0.5)]);    
    replacements{end+1} = {'C1_outer_ind',             commavec(outerStarts-1)};
    replacements{end+1} = {'C1_inner_ind',             commavec(innerInd-1)};        
end



for i = 1:length(replacements)
    source = replacements{i}{1};
    target = replacements{i}{2};
    if isa(target,'double') || isa(target,'logical')
        target = num2str(target);
    end
    
    code_template_str = strrep(code_template_str,['$(' source ')'], target);
end

template.code_str = code_template_str;

if opt.write_to_file
    if ~exist('./solvers_cpp','dir')
        mkdir('solvers_cpp')
    end

    fname = ['solvers_cpp/solver_' solv_name '.cpp'];
    fprintf('Writing to "%s"...',fname);
    fid = fopen(fname,'w');
    fprintf(fid,'%s\n',template.code_str);
    fclose(fid);
    fprintf(' OK\n');
    
    fprintf('Building "%s"...\n',fname);
    if opt.cg_compile_mex
        
        if strcmp(opt.eigen_solver,'sturm') || strcmp(opt.eigen_solver,'sturm_dani')
            mex(['-I"' opt.cg_eigen_dir '"'],['-I"' template_dir '"'],'-O', '-R2018a', fname,'-outdir','solvers_cpp')
        else
            mex(['-I"' opt.cg_eigen_dir '"'],'-O', '-R2018a', fname,'-outdir','solvers_cpp')
        end
    end
    
    template.solv_fun = str2func(['solver_' solv_name]);
else
    template.solv_fun = [];
end

function str = commavec(v)
str = sprintf('%d,',v);
str = str(1:end-1);

function str = print_debug_comments(solv,template,opt)
% Debug
tmp = 1;
if ~isempty(opt.saturate_mon)
    tmp = template.satmon;
end
str = sprintf('\n%s%s\n','// Action = ', sym2char(sym(solv.actions).'));
str = [str sprintf('%s%s\n','// Quotient ring basis (V) = ', sym2char(sym(solv.basis) ./tmp))];
available_mons = [template.reducible template.basis];
str = [str sprintf('%s%s\n','// Available monomials (RR*V) = ', sym2char(sym(available_mons) ./tmp))];



