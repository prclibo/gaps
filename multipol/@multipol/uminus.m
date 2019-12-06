function p1 = uminus(p1)
% MULTIPOL/UMINUS operator
% -p1 negates the coefficients of the polynomials
% p2 = uminus(p1);
% p2 = -p1;

for i=1:numel(p1)
	p1(i).coeffs = -p1(i).coeffs;
end
