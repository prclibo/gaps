function n = nterms(p)
%returns the number of terms in a multipol object
n = zeros(size(p));
for i=1:numel(p)
	n(i) = numel(p(i).coeffs);
end