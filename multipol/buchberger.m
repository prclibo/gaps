function G = buchberger(p,order,maxiter,tol_min,tol_max)
% Simple implementation of Buchberger's algorithm
if nargin<4
	tol_min = 100*eps;
	tol_max = inf;
end
G = eqsize(p(:));
zero = multipol;
iter = 0;
while true
	iter = iter+1;
	G2 = G;
	for i=1:numel(G2)-1
		for j=i+1:numel(G2)
			fprintf('.');
			S = spoly(G2(i),G2(j),order,tol_min)/G;
			S = squeeze(S,tol_min);
			if S~=zero
				G = [G; S]; %#ok
			end
		end
		fprintf(' ');
	end
	
	G = unique(G);
	[~,ia] = intersect(G,-G); % Also remove polynomials duplicate up to sign
	if numel(ia)>0
		G(ia(1:end/2)) = [];
	end
	
	rem = false(numel(G),1);
	for i=1:numel(G)
		if all(abs(coeffs(G(i)))>tol_max) || all(abs(coeffs(G(i)))<tol_min), rem(i) = true; end
	end
	G(rem) = [];
	
	fprintf('\n%u: %u polynomials\n',iter,numel(G));
	if iter==maxiter, break; end
	if numel(G)==numel(G2)
		break;
	end
end