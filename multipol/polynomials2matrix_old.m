function [M occurringMonomials] = polynomials2matrix_old(p,order)
% [M occurringMonomials] = polynomials2matrix(polynomials)
% Creates the matrix of coefficients obtained by stacking
% the polynomials in 'polynomials' on top of each other
if nargin<2
	order = 'grevlex';
end

p = eqsize(p);

if(iscell(p))
    [M occurringMonomials] = polynomials2matrix_cell(p);
    return;
end

% assemble occurring monomials
occurringMonomials = inf(1,nvars(p(1)));
for i = 1 : length(p)
  occurringMonomials = union(monomials(p(i))', occurringMonomials, 'rows');
end
occurringMonomials = occurringMonomials(1:end-1, :)';
occurringMonomials = multipol(ones(1, size(occurringMonomials, 2)),...
    occurringMonomials);
occurringMonomials = sort(occurringMonomials,order);

% determine matrix size
rows = length(p);
cols = nterms(occurringMonomials);

% create matrix
M = zeros(rows,cols);
for i = 1 : length(p)
    c = coeffsof_p2m(occurringMonomials, p(i));
	M(i,:) = c;
end

% create monomial vector
occurringMonomials = mons2vec(occurringMonomials);

% remove zero columns
ind = sum(abs(M),1)~=0;
M = M(:, ind);