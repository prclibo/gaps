function [M occurringMonomials] = polynomials2matrix_cell(polynomials)
% N.B! This is an old method kept for compatibility. Use polynomials2matrix
% instead.
%
% [M occurringMonomials] = polynomials2matrix_cell(polynomials)
% Creates the matrix of coefficients obtained by stacking
% the polynomials in 'polynomials' on top of each other

if(length(polynomials)==1)
    tmp = polynomials;
    polynomials = cell(1,1);
    polynomials{1} = tmp;
end

%assemble occurring monomials
occurringMonomials = cell(size(polynomials));
for i = 1 : length(polynomials)
  occurringMonomials{i} = monomials(polynomials{i});
end
occurringMonomials = cell2mat(occurringMonomials);
occurringMonomials = sort(multipol(ones(1,size(occurringMonomials,2)),occurringMonomials));
occurringMonomials = unitycoeffs(occurringMonomials);
occurringMonomials = sorttlex(occurringMonomials);

%determine matrix size
rows = length(polynomials);
cols = nterms(occurringMonomials);

% create matrix
M = zeros(rows,cols);
for i = 1 : length(polynomials)
    c = coeffsof_p2m(occurringMonomials, polynomials{i});
    M(i,:) = c;
end

% remove zero columns
ind = find(sum(abs(M),1));
M = M(:, ind);





