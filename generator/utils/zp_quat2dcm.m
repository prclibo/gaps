function R = zp_quat2dcm(q, p)

assert(logical(mod(sum(q.^2), p) == 1));
R = sym_quat2dcm(q);
R = mod(R, p);

