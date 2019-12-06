function [nonPermissible permissible] = getPermissible_old(mon, variable)
% N.B. This function is kept for compability. Use getPermissible instead!
%
%[nonPermissible permissible] = getPermissible_old(mon, variable) splits the
% set of monomials in mon into two sets based on multiplication with the
% variable indexed in "variable". if "variable" is negative, division
% instead of multiplication is considered.

mm = monomials(mon);
n = size(mm, 1);
var = zeros(n, 1);
var(abs(variable)) = sign(variable);
permissible = [];
nonPermissible = [];

for k = 1 : size(mm, 2)
    m = mm(:, k);
    m = m + var;
    comp = mm ~= repmat(m, 1, size(mm, 2));
    if(sum(sum(comp) == 0) > 0)
        permissible = [permissible k];
    else
        nonPermissible = [nonPermissible k];
    end
end