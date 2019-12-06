function m = monomials(p1);
% MULTIPOL/MONOMIALS operator
% Returns the monomials of the polynom p1
% m = monomials(p1);

m={p1.monomials};
if numel(m) == 1, m = m{1}; end
