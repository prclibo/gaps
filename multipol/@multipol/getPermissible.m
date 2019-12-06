function [nonPermissible permissible] = getPermissible(mon, variable)
% [nonPermissible permissible] = getPermissible(mon, variable) splits the
% set of monomials in mon into two sets based on multiplication with the
% variable indexed in "variable". If "variable" is negative, division
% instead of multiplication is considered. If "variable" set to 'all' then
% "permissible" will be an intersection of all permissible sets for
% different variables.

permissible = 1 : length(mon);
nonPermissible = [];
p = nvars(mon(1));

if(strcmp(variable, 'all'))
    for k = 1 : p
        [non_perm_k perm_k] = getPermissible(mon, k);
        permissible = intersect(permissible, perm_k);
        nonPermissible = union(nonPermissible, non_perm_k);
    end
    return;
end

if(length(mon) == 1)
    [nonPermissible permissible] = getPermissible_old(mon, variable)
    return;
end

mm = monvec2matrix(mon);
n = size(mm, 1);
var = zeros(n, 1);
var(abs(variable)) = sign(variable);
permissible = [];
nonPermissible = [];

for k = 1 : size(mm, 2)
    m = mm(:, k);
    m = m + var;
    comp = mm ~= repmat(m, 1, size(mm, 2));
    if(sum(sum(comp, 1) == 0) > 0)
        permissible = [permissible k];
    else
        nonPermissible = [nonPermissible k];
    end
end