function p_out = divide_by_1_plus_x2(p);
% p_out = divide_by_1_plus_x2(p);
% polynomial division for the special case of division by (x^2 + 1)
% where x is the first variable in the variable ordering

nVariables = size(p.monomials, 1);

%define the x
x = multipol(1,[1; zeros(nVariables - 1, 1)]);
x2 = x*x;

p = sort(p);
p_out = multipol(0,zeros(nVariables, 1));
monoms = fliplr(p.monomials);
maxdegree = monoms(1,end);

k=1;
%iterate through the terms in p
while(true)
    %extract coefficients and monomials from the remainder polynomial
    coeffs = p.coeffs;
    coeffs = fliplr(p.coeffs);
    monoms = p.monomials;
    monoms = fliplr(p.monomials);
    
    %terminate if there are no more terms in the remainder polynomial
    if(k>length(coeffs)); break; end;
    
    %get the monomial to eliminate
    currentMonomial = multipol(coeffs(:,k), monoms(:,k));
    
    %terminate if we wont be able to eliminate
    if(monoms(1,k) > maxdegree - 2); break; end;
    
    %build quotient polynomial
    p_out = p_out + currentMonomial;
    
    %substract from remainder polynomial
    currentMonomial = currentMonomial * x * x;
    p = p - currentMonomial;
    
    k=k+1;
end

%remove zero coeffs
p_out = squeeze(p_out);