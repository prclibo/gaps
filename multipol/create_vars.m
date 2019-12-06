function x = create_vars(n,perm)
% Return n multipol variables
% To achieve different monomial orderings, this
% is a good place to permute the underlying variables

if nargin<2
	perm = 1:n;
end

for i=n:-1:1
	z = zeros(n,1);
	z(perm(i)) = 1;
	x(i,1) = multipol(1,z);
end