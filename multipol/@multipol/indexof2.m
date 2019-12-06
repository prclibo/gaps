function [ind2 ind] = indexof2(monoms, polynomial)
%ind = indexof2(monoms, polynomial)
%returns a vector of indices of the monomials in "monoms" indicating the
%position in "polynomial". Both monoms and polynomial are multipol objects.

%keyboard;
n = nterms(polynomial);
c = 1 : n;
p = polynomial;
p.coeffs = c;
ind = coeffsof(monoms, p);
[blubb,sorti]=sort(ind);
nz = find(blubb>0);
ind2 = sorti(nz);



