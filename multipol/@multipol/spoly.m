function S = spoly(f,g,order,tol)
% Compute the S-polynomial of f and g
if nargin<3
	tol = 100*eps;
end

f = sort(f,order);
g = sort(g,order);
LMf = f.monomials(:,1);
LMg = g.monomials(:,1);
LCM = multipol(1,max(LMf,LMg));

LTf = multipol(1,LMf)*f.coeffs(1);
LTg = multipol(1,LMg)*g.coeffs(1);

f = multipol(f.coeffs(2:end),f.monomials(:,2:end)); % Remove first term to ensure it cancels
g = multipol(g.coeffs(2:end),g.monomials(:,2:end));

S = LCM * (f./LTf - g./LTg);
S = squeeze(S,tol);