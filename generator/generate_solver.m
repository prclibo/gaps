function [solv,opt] = generate_solver(prob_fn, opt)

prob = prob_fn();

solv.prob = prob;
solv.name = class(prob);
solv.opt = opt;
solv.n_vars = numel(prob.unk_vars);
solv.vars = multipol(prob.unk_vars, 'vars', prob.unk_vars);
solv.eqs_zp = cell(1, opt.integer_expansions);
solv.unk_zp = cell(1, opt.integer_expansions);
fprintf('Sampling Zp instances\n'); tic;
for i = 1:numel(opt.integer_expansions)
    solv.eqs_zp{i} = prob.rand_eq_zp(opt.prime_field);
end
toc; fprintf('-- Sampled.\n'); tic;

solv.n_eqs = numel(prob.eqs);



%%
fprintf('-------------------------------------------------\n');
fprintf('generate_solver'); % (%s, %d equations in %d variables)\n',solv.name,solv.n_eqs,solv.n_vars);


%% Check for symmetries in the problem

if opt.find_sym
    % TODO Transplant
    fprintf('Checking for symmetry. max_p = %d\n',opt.sym_max_p); tic;
    
    [opt.sym_cc,opt.sym_pp] = find_symmetries(solv.eqs_zp{1},opt.sym_max_p);
    
    toc; fprintf('-- Found %d symmetries.\n',length(opt.sym_pp));
    for k = 1:length(opt.sym_pp)
        fprintf('\tc = [ ');
        fprintf('%d ',opt.sym_cc(:,k));
        fprintf('],\tp = %d\n',opt.sym_pp(k));
    end
end



%% Compute monomial basis for quotient space
if isempty(opt.custom_basis)
    solv.basis = find_monomial_basis(solv.eqs_zp{1}, solv,opt,solv.name);
else
    solv.basis = opt.custom_basis;
    fprintf('Using custom basis for quotient space.\n');
end
solv.n_sol = length(solv.basis);

fprintf('Problem has (at most) %d solutions.\n', solv.n_sol);
fprintf('Quotient ring basis (degree = %d)\n', length(solv.basis));
fprintf('\tB = [ %s ]\n', strjoin(string(solv.basis)));

%% Remove zero solution
if opt.remove_zero_sol && ...
        all(evaluate(solv.eqs_zp{1},zeros(solv.n_vars,1)) == 0)
    
    fprintf('Removing zero solution.\n');
    
    % try to remove constant from quotient basis
    % (it must be independent from the rest since zero is in variety)
    i = find(solv.basis == multipol(1));
    if ~isempty(i)
        solv.basis(i) = [];
    end
end

%% Select action monomial

if isempty(opt.actmon)
    solv.actions = select_actions(solv,opt);
else
    solv.actions = opt.actmon;
end
fprintf('Action monomial = [ %s ]\n', strjoin(string(solv.actions)));

if opt.use_sym
    % make sure action monomial is invariant
    for i = 1:length(solv.actions)
        for k = 1:length(opt.sym_pp)
            if mod(opt.sym_cc(:,k)'*monomials(solv.actions(i)),opt.sym_pp(k))~=0
                error('Action monomial must be invariant under the symmetry!');
            end
        end
    end
end

%% Select monomials to reduce

tic
solv.reducible = [];
for k = 1:length(solv.actions)
    solv.reducible = [solv.reducible solv.actions(k)*solv.basis];
end
if opt.force_vars_in_reducibles
    solv.reducible = unique([solv.reducible, solv.vars]);
end
if ~isempty(opt.extra_reducible)
    solv.reducible = unique([solv.reducible opt.extra_reducible]);
end
toc
% remove basis elements
solv.reducible = mono_setdiff(solv.reducible, solv.basis);

fprintf('Monomials to reduce: (%d monomials)\n',length(solv.reducible));
fprintf('\tR = [ %s ]\n', strjoin(string(solv.reducible), ', '));

%% Generate monomial expansion

solv.templates = initialize_templates(solv,opt);

fprintf('Finding elimination templates ... ');
As = cell(1,length(solv.eqs_zp));
for k = 1:length(solv.eqs_zp)
    id = [solv.name '.' num2str(k)];
    [As{k},opt] = find_template(solv.eqs_zp{k},solv,opt,id);
end
fprintf('OK\n');

% merge expansions
A = cell(size(As{1}));
for k = 1:numel(A)
    tmp = cellfun(@(x) x{k},As,'UniformOutput',0);
    A{k} = cat(2,tmp{:});
    A{k} = -unique(-A{k}','rows')';
end


solv = build_templates(A,solv,opt);


%%

% Select which templates to construct solvers from
if opt.build_all_templates
    solv.target_template = true(length(solv.templates),1);
else
    % only build best template
    [~,template_ind] = sortrows([cellfun(@(x) x.sz, solv.templates);cellfun(@(x) length(x.reducible), solv.templates)]');    
    solv.target_template = false(length(solv.templates),1);
    solv.target_template(template_ind(1)) = 1;
end

% Print results
fprintf('Found %d elimination templates.\n',length(solv.templates));
desc_len = max(cellfun(@(x) length(x.desc),solv.templates));
fmt = ['%s %-' num2str(desc_len+1) 's - %- 4d rows, %- 2d basis, %- 2d reducible, %s\n'];
extra = '';
for k = 1:length(solv.templates)
    if solv.target_template(k), build_mark = '*'; else, build_mark = ' '; end
    if ~isempty(opt.saturate_mon), extra = sprintf('N = %d',solv.templates{k}.saturate_degree);end

    fprintf(fmt,build_mark,solv.templates{k}.desc,solv.templates{k}.sz,...
        length(solv.templates{k}.basis),length(solv.templates{k}.reducible),extra);
end

if opt.stop_after_template
    return;
end

%% Compute coefficients
fprintf('Extracting template coefficients... ');
solv.coefficients = compute_coefficients(prob, solv,opt);
fprintf('OK (%d found)\n',length(solv.coefficients.coeff_eqs));



for tk = find(solv.target_template)'
    fprintf('Building solver (%s_%s)\n',solv.name,solv.templates{tk}.desc);   
       
    tic
    solv.templates{tk} = finalize_template(solv,solv.templates{tk}, opt);
    tt = toc;

    fprintf('Elimination template size = [%d,%d], nnz = %d\n',solv.templates{tk}.C_sz,length(solv.templates{tk}.C_ind));

    if opt.remove_extra_columns
        solv.templates{tk} = remove_extra_columns(solv.templates{tk}, solv);
    end


    if opt.build_all_templates
        name = [solv.name '_' solv.templates{tk}.desc];
    else
        name = solv.name;
    end
    
    if strcmp(opt.cg_language,'matlab')
        solv.templates{tk} = generate_code(name,solv,solv.templates{tk},opt);
    elseif strcmp(opt.cg_language,'cpp_eigen') || strcmp(opt.cg_language,'cpp')
        solv.templates{tk} = generate_code_cpp(name,solv,solv.templates{tk},opt);
    else
        error('Unknown output language');
    end
end

return
