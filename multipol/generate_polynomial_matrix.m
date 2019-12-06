function [M, resultMonomials] = generate_polynomial_matrix(polynomials, productmonomials)
% [M,resultMonomials] = generate_polynomial_matrix(polynomials,
% productMonomials)
%
% generates a matrix by multiplying all the polynomials with the product
% monomials and then stacking the coefficients.

newPolys = {};
p=1;
for k = 1:length(productmonomials);
    monoms = monomials(productmonomials{k});
    for m = 1 : size(monoms, 2)
        mon = multipol(1, monoms(:,m));
        newPolys{p} = polynomials{k}*mon;
        p = p + 1;
    end
end

[M resultMonomials] = polynomials2matrix(newPolys);