function p1 = remove_monoms(p1,indeces);
% MULTIPOL/REMOVE_MONOMS
% Removes the monom on the places given in indeces

p1.monomials(:,indeces) = [];
p1.coeffs(:,indeces) = [];

