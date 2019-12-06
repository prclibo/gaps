function newCoeffs = polynomial2vector(polynomial)
%c = polynomial2vector(polynomial)
%generates a vector with the coefficients of 'polynomial'.
%unlike polynomial.coeffs this returns the coefficients for
% *all* monomials, i.e. possibly lots of zeros.

%Note SiBu 2012-02-17
%Ok should work now. fixdit

polynomial = sorttlex(polynomial);
m = polynomial.monomials;
c = polynomial.coeffs;

%define degrees of monomials to create
maxdegree = sum(m(:,1));
maxleading = max(m(:,1));
nVariables = size(m,1);
degrees = ones(1, nVariables) * maxdegree;

%create monomials
newMonomials = create_upto_degree(degrees, maxdegree); 
newMonomials = newMonomials + multipol(1,0); 
newMonomials = sorttlex(newMonomials);
newMonomials = newMonomials.monomials;


%set the coefficients

newCoeffs = zeros(1, size(newMonomials, 2));
for i = 1 : size(m, 2)
    mon = m(:, i);
    mon = mon * ones(1, size(newMonomials, 2));
    logical = double(mon==newMonomials);
    ind = find(prod(logical, 1));
    newCoeffs(ind) = c(i);
end


%reorder in ascending order
%newCoeffs = fliplr(newCoeffs);

end
