function p = setcoeff(p, monomial, value)
% p = setcoeff(p, monomial, value) sets the coefficient of monomial to
% value.
%INPUT
%p          Multipol object   polynomail to be changed coeffecient
%monomial   multipol object
%value      scalar New coefficient value

m = monomials(p);
c = coeffs(p);
mon = repmat(monomials(monomial), 1, length(c));
ind = find(prod((mon == m)*1));
c(ind) = value;
p = multipol(c, m);