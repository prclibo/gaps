function ind = indexof(monoms, polynomial)
%ind = indicesof(monoms, polynomial)
%returns a vector of indices of the monomials in "monoms" indicating the
%position in "polynomial". Both monoms and polynomial are multipol objects.

keyboard;
n = size(polynomial);
c = 1 : n;
p = polynomial;
p.coeffs = c;
ind = coeffsof(monoms, p);

