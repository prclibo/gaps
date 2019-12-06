function p = squeeze(p,tol)
% MULTIPOL/SQUEEZE operator
% Squeezes out terms with zero coefficient
% p2 = squeeze(p1);

if nargin<2
	tol = 100*eps;
end

% if isempty(p.monomials)
% 	p.monomials = [];
% 	if isempty(p.coeffs) || abs(p.coeffs)<tol
% 		p.coeffs = 0;
% 	end
% 	return
% end

for i=1:numel(p)
	
	zero = find(abs(p(i).coeffs)<tol);
	p(i).coeffs(zero) = [];
	p(i).monomials(:,zero) = [];
	
	if isempty(p(i).coeffs)
		p(i).coeffs = 0;
		p(i).monomials = zeros(size(p(i).monomials,1),1);
	end
	
end