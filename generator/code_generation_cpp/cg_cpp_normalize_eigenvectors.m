function str = cg_cpp_normalize_eigenvectors(solv,template,opt)
str = [];
if isempty(template.norm_scheme)
    str = [str '#error TODO: Normalize eigenvectors.\n'];
elseif strcmp(template.norm_scheme.method,'direct')
    % We have 1 in basis.
    str = [str sprintf('V = (V / %s.replicate(%d, 1)).eval();\n\n',...
        cg_cpp_parse_source(template.norm_scheme.sources),...
        length(template.basis))];
elseif strcmp(template.norm_scheme.method,'eigen')
    mm = cg_cpp_parse_scheme(template.norm_scheme);
    k = length(template.norm_scheme.sources);        
    expr = ['D.transpose() / ' mm];
    
    if k == 2
        str = [str sprintf('ArrayXXcd scale = (%s).sqrt();\n',expr)];
    elseif k > 2
        str = [str sprintf('ArrayXXcd scale = (%s).pow(1/%d);\n',expr,k)];
    else
        str = [str sprintf('ArrayXXcd scale = %s;\n',expr)];
    end    
    str = [str sprintf('V = (V * scale.replicate(%d, 1)).eval();\n\n',...
        length(template.basis))];
elseif strcmp(template.norm_scheme.method,'sqrt_eigen')
    mm = cg_cpp_parse_scheme(template.norm_scheme);
    k = length(template.norm_scheme.sources);        
    expr = ['(D.transpose().sqrt() / ' mm ')'];
    if k == 2
        str = [str sprintf('ArrayXXcd scale = %s.sqrt();\n',expr)];
    elseif k > 2
        str = [str sprintf('ArrayXXcd scale = %s.pow(1/%d);\n',expr,k)];
    else
        str = [str sptrinf('ArrayXXcd scale = %s;\n',expr)];
    end
    str = [str sprintf('V = (V * scale.replicate(%d, 1)).eval();\n\n',...
        length(template.basis))];
end

