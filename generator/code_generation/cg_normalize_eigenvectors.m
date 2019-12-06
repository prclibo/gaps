function str = cg_normalize_eigenvectors(solv,template,opt)
str = [];
if isempty(template.norm_scheme)
    str = [str sprintf('error(''TODO: Normalize eigenvectors.'');\n')];
elseif strcmp(template.norm_scheme.method,'direct')
    % we have 1 in basis
    str = [str sprintf('V = V ./ (ones(size(V,1),1)*%s);\n',cg_parse_source(template.norm_scheme.sources))];
elseif strcmp(template.norm_scheme.method,'eigen')   
    mm = cg_parse_scheme(template.norm_scheme);
    k = length(template.norm_scheme.sources);        
    expr = ['(diag(D).'' ./ (' mm '))'];
    
    if k == 2        
        str = [str sprintf('scale = sqrt%s;\n',expr)];
    elseif k > 2
        str = [str sprintf('scale = %s.^(1/%d);\n',expr,k)];
    else
        str = [str sprintf('scale = %s;\n',expr)];
    end    
    str = [str sprintf('V = V .* (ones(size(V,1),1)*scale);\n')];
elseif strcmp(template.norm_scheme.method,'sqrt_eigen')   
    mm = cg_parse_scheme(template.norm_scheme);
    k = length(template.norm_scheme.sources);        
    expr = ['(sqrt(diag(D)).'' ./ (' mm '))'];
    if k == 2        
        str = [str sprintf('scale = sqrt%s;\n',expr)];
    elseif k > 2
        str = [str sprintf('scale = %s.^(1/%d);\n',expr,k)];
    else
        str = [str sptrinf('scale = %s;\n',expr)];
    end    
    str = [str sprintf('V = V .* (ones(size(V,1),1)*scale);\n')];
end

