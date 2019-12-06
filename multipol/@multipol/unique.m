function [u ia ic] = unique(eq)
% MULTIPOL/UNIQUE
eq = eqsize(eq);

c1 = cell(size(eq));
for i=1:numel(eq)
	c1{i} = hash(eq(i));
end

[~,ia,ic] = unique(c1,'stable');
ia = ia(:)';
ic = ic(:)';
u = eq(ia);