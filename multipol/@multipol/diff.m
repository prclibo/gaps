function J = diff(eq,var)
% MULTIPOL/DIFF Jacobian of polynomials in eq or derivative wrt variables
% in "var"

if nargin<2
	var = 1:nvars(eq(1));
end

for v=numel(var):-1:1
	d = eq(:);
	for k=numel(eq):-1:1
		d(k).coeffs = d(k).coeffs.*d(k).monomials(var(v),:);
		d(k).monomials(var(v),:) = d(k).monomials(var(v),:)-1;
		d(k) = squeeze(d(k));
	end
	J(:,v) = d;
end