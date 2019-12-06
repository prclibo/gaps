function a_zp = zp_eval_rational(a, p)

[num, den] = numden(a);
num = mod(num, p);
den = mod(den, p);
inv_den = zp_inv(den, p);

a_zp = mod(num .* inv_den, p);