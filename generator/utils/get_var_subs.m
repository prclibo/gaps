function var_subs = get_var_subs(symb, var)

dims = nnz(size(var) > 1);

if dims >= 1
    sz = size(var);
    nm = [symb, repmat('%d', 1, dims)];
    var_subs = struct('var', num2cell(sym(nm, size(var))),...
        'subs', num2cell(var));
    var_subs = var_subs(:);
else
    var_subs = struct('var', num2cell(sym(symb)),...
        'subs', num2cell(var));
end