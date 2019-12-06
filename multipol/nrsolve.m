function [x0 res] = nrsolve(eq,x0,maxiter)

eq = eqsize(eq(:));
n = nvars(eq(1));

if nargin<3
	maxiter = 50;
end

if nargin<2
	x0 = randn(n,1)*10 + 1i*randn(n,1)*10;
end

J = diff(eq);
res = zeros(1,size(x0,2));
for k=1:size(x0,2);
	
	x = x0(:,k);

	for iter=1:maxiter

		fx = evaluate(eq,x); if iter==1, res0 = norm(fx); end
	% 	norm(fx)
		
		if norm(fx)<1e-14, break; end
		Jx = evaluate(J,x);

		nx = x - Jx\fx;
% 		nx = x - pinv(Jx)*fx;
% 		nx = x - (Jx'*Jx+eye(n)*1e-6)\(Jx'*fx);
	% 	if norm(nx-x)<1e-9, break; end
		x = nx;

	end
% 	log10(svd(Jx))'
% 	Jx
	% iter
	res(k) = norm(evaluate(eq,x));
	fprintf('Iter.: %2u, residual: %-7g --> %g\n',iter,res0,res(k));
% 	cond(Jx)
% 	if rank(Jx)<min(size(Jx)), keyboard; end
	x0(:,k) = x;
	
end


end