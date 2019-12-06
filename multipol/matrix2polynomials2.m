function polys = matrix2polynomials2(matrix, mons)
% polys = matrix2polynomials(matrix, mons) creates a cell array of multipol
% polynomials from the coefficients in matrix using the monomials in the
% matrix mons.

for k = 1 : size(matrix, 1)
    polys(k,1) = multipol(matrix(k, :), mons);
end
