function ret = zp_normr(vec, p)

sqr_norm = mod(sum(powermod(vec, 2, p)), p);
len = zp_sqrt(sqr_norm, p);
assert(~isempty(len));

inv_len = zp_inv(len(1), p);
ret = mod(vec * inv_len, p);