function [p, ia] = mono_setdiff(eq1, eq2)

mdeg1 = monomials(eq1);
mdeg2 = monomials(eq2);

if ~iscell(mdeg1), mdeg1 = {mdeg1}; end
if ~iscell(mdeg2), mdeg2 = {mdeg2}; end

nterm1 = cellfun(@(x) size(x, 2), mdeg1);
nterm2 = cellfun(@(x) size(x, 2), mdeg2);
assert(all(nterm1 == 1) && all(nterm2 == 1), 'Only works for mono but not poly');

mdeg1 = cat(2, mdeg1{:});
mdeg2 = cat(2, mdeg2{:});

[~, ia] = setdiff(mdeg1.', mdeg2.', 'rows', 'stable');
ia = ia(:).';
p = eq1(ia);

end

