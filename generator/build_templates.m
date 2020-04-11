function [ solv ] = build_templates( A, solv, opt )

for k = 1:length(solv.templates)
    solv.templates{k} = build_template(solv.templates{k},A,opt);
end

function template = build_template(template,A,opt)

    % which reducibles do we need
    basis = template.basis;
    reducible = template.reducible;

    A = A(:,template.reducible_ind);

    % Take care of saturation
    if ~isempty(opt.saturate_mon)
        degs = opt.saturate_degree(template.reducible_ind);        
        template.saturate_degree = max(degs);
        template.satmon = opt.saturate_mon^template.saturate_degree;

        basis = basis * template.satmon;
        reducible = reducible * template.satmon;

        % if some reducible was a lower degree, multiply coefficients
        for k = 1:length(degs)
            if degs(k) < template.saturate_degree
                tmp = monomials(opt.saturate_mon^(template.saturate_degree-degs(k)));
                for i = 1:size(A,1)
                    A{i,k} = A{i,k} + tmp * ones(1,size(A{i,k},2));
                end
            end
        end
    end
    
    % template.monomials{i}: all monos need to multiple on eq_i
    template.monomials = cell(size(A,1),1);    
    for i = 1:size(A,1)
        Ai = A(i,:);        
        template.monomials{i} = unique(cat(2,Ai{:})','rows')';
        if isempty(template.monomials{i})
            template.monomials{i} = zeros(nvars(basis(1)), 1);
        end
    end
    template.sz = sum(cellfun(@(x) size(x,2),template.monomials));  
   



