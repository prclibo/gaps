function unit = zp_rand_unit(d, p)

while true
    vec = sym(randi(p - 1, d, 1));
    sqr_norm = mod(sum(powermod(vec, 2, p)), p);
    len = zp_sqrt(sqr_norm, p);
    if ~isempty(len)
        inv_len = zp_inv(len(1), p);
        unit = mod(vec * inv_len, p);
        return
    end
end