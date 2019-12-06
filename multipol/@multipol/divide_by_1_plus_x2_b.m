function [quotient remainder] = divide_by_1_plus_x2(p);
% [quotient remainder] = divide_by_1_plus_x2(p);
% polynomial division for the special case of division by (x^2 + 1)
% where x is the first variable in the variable ordering.

mons = monomials(p);
maxdegree = max(mons(1,:));
remainder = p;
dim = size(mons, 1);
quotient = multipol(0, [zeros(dim,1)]);
x = multipol(1,[1;zeros(dim-1,1)]);
x2 = multipol(1,[2;zeros(dim-1,1)]);

while(maxdegree >= 2)
  %get the coefficient of the leading term a_max
  a_max = coeffpol(remainder, x^maxdegree);
  
  %add a_max*x^(maxdegree - 2) to the quotient polynomial
  quotient = quotient + a_max*x^(maxdegree - 2);
  
  %subtract a_max*(x^2 + 1)*x^(maxdegree - 2) from the remainder
  %polynomial
  remainder = remainder - a_max*(x^2 + 1)*x^(maxdegree - 2);
  
  %get the new max degree
  mons = monomials(remainder);
  maxdegree = max(mons(1,:));
end

%remove zero coeffs
quotient = squeeze(quotient);
remainder = squeeze(remainder);
