function template = generate_code(solv_name, solv,template,opt)

if ~isfield(template,'C_sz') || isempty(template.C_sz)
    fprintf('Template indicies are missing. Generating...\n');
    template = finalize_template(solv,template);
end

addpath([fileparts(mfilename('fullpath')) filesep 'code_generation'])

prob = solv.prob;
in_par_str = strjoin(fieldnames(prob.in_subs), ', ');
out_par_str = strjoin(fieldnames(prob.out_subs), ', ');
str = sprintf('function [%s] = solver_%s(%s)\n', out_par_str, solv_name, in_par_str);

% Construct template C
str = [str sprintf('[C0,C1] = setup_elimination_template(%s);\n', in_par_str)];

% Construct action matrix AM
str = [str cg_setup_action_matrix(solv,template,opt)];

% Compute eigenvectors
if ~opt.generalized_eigenvalue_solver
    if opt.sparse_template
        str = [str sprintf('[V,D] = eig(full(AM));\n')];    
    else
        str = [str sprintf('[V,D] = eig(AM);\n')];
    end
else
    str = [str sprintf('[V,D] = eig(A,B);\n')];
end

% Fix the scale on the eigenvectors
str = [str cg_normalize_eigenvectors(solv,template,opt)];

% Extract solution from eigenvectors
str = [str cg_extract_solutions(solv,template,opt)];

str = [str print_debug_comments(solv,template,opt)];

% Add compute coeff function
str = [str cg_compute_coeff(solv,template,opt)];

% Add template function
str = [str cg_setup_template(solv,template,opt)];

template.code_str = str;

if opt.write_to_file
    if ~exist('./solvers','dir')
        mkdir('solvers')
    end

    fname = ['solvers/solver_' solv_name '.m'];
    fprintf('Writing to "%s"...',fname);
    fid = fopen(fname,'w');
    fprintf(fid,'%s\n',template.code_str);
    fclose(fid);
    fprintf(' OK\n');
    template.solv_fun = str2func(['solver_' solv_name]);
else
    template.solv_fun = [];
end


function str = print_debug_comments(solv,template,opt)
% Debug
tmp = 1;
if ~isempty(opt.saturate_mon)
    tmp = template.satmon;
end
str = sprintf('\n%s%s\n','% Action = ', sym2char(sym(solv.actions).'));
str = [str sprintf('%s%s\n','% Quotient ring basis (V) = ', sym2char(sym(solv.basis) ./tmp))];
available_mons = [template.reducible template.basis];
str = [str sprintf('%s%s\n','% Available monomials (RR*V) = ', sym2char(sym(available_mons) ./tmp))];



