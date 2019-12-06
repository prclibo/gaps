function mon = monvec(p)
% mon = monvec(p) returns the vector of monomials occurring in p.
m = monomials(p);
for k = 1 : size(m, 2);
    mon(k) = multipol(1, m(:, k));
end
mon = mon(:);