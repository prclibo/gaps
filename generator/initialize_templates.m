function [ templates ] = initialize_templates( solv, opt )

templates = cell(1,length(solv.actions));

% which reducibles should be included in all templates?
all_inc = false(1,length(solv.reducible));
if opt.force_vars_in_reducibles
    red_vars = mono_setdiff(solv.vars, solv.basis);
    all_inc(find_mon_indices(solv.reducible, red_vars)) = 1;
end
if ~isempty(opt.extra_reducible)
    red_vars = setdiff(opt.extra_reducible, solv.basis);
    all_inc(find_mon_indices(solv.reducible, red_vars)) = 1;
end

if opt.use_sym && ~isempty(opt.sym_pp)
    templates = initialize_templates_symmetry(solv, opt, all_inc);
else
    templates = cell(1,length(solv.actions));
    % for each action, figure out the template
    for k = 1:length(solv.actions)
        action = solv.actions(k);
        templates{k} = initialize_template(action, solv.basis, solv.reducible, all_inc,opt);
    end
end
end


function templates = initialize_templates_symmetry(solv, opt, all_inc)
  
    B = monvec2matrix(solv.basis);

    tmp = opt.sym_cc'*B;
    for k = 1:size(tmp,1)
        tmp(k,:) = mod(tmp(k,:),opt.sym_pp(k));
    end

    if isempty(opt.sym_rem)
        % generate possible remainders
        rems = unique(tmp','rows')';
    else
        rems = opt.sym_rem(:);
    end
    
    templates = {};
    % for each symmetry we split into the remainders
    for k = 1:size(rems,2)

        % select basis and reducible for this remainder
        basis_ind = sum(abs(tmp - rems(:,k) * ones(1,size(tmp,2))),1)==0;
        
        basis_k = solv.basis(basis_ind);

        for i = 1:size(solv.actions)
            templates{end+1} = initialize_template(solv.actions(i),basis_k,solv.reducible,all_inc,opt);
            templates{end}.sym_rem = rems(:,k);
            templates{end}.desc = sprintf('%s_rem%s',templates{end}.desc,sprintf('_%d',rems(:,k)));
        end
    end
end

function template = initialize_template(action, basis, reducible, all_inc, opt)

    % which reducibles do we need
    inc = all_inc;
    inc(find_mon_indices(reducible,mono_setdiff(action*basis,basis))) = 1;
    



    template.reducible_ind = inc;    
    template.monomials = [];
    
    template.sz = nan;  
    template.action = action;
    template.basis = basis;
    template.reducible = reducible(inc);
    template.desc = strrep(strrep(sprintf('action_%s',string(action)),'^','__'),'*','');
    template.extraction_scheme = [];
    template.norm_scheme = [];

end

