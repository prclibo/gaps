function n = nvars(p)
%returns the number of variables in a multipol object
n = zeros(size(p));
for i=1:numel(p)
	n(i) = size(p(i).monomials,1);
end