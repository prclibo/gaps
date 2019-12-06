function R = sym_quat2dcm(q)


s = q(1);
u1 = q(2);
u2 = q(3);
u3 = q(4);
u = [u1; u2; u3];
R = 2 * (u * transpose(u) - s * skew(u)) + (s * s - transpose(u) * u) * eye(3);
