function [ ainv ] = zp_inv( a, p )

ainv = arrayfun(@(x)zp_inv_one(x, p), a);

end
function [ainv] = zp_inv_one(a, p)
% https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm

t = sym(0); new_t = sym(1);
r = p; new_r = a;
while logical(new_r ~= 0)
    q = floor(r/new_r);
    [t, new_t] = deal(new_t, t - q * new_t);
    [r, new_r] = deal(new_r, r - q * new_r);
end
if logical(t < 0)
    t = t + p;
end

ainv = t;
end