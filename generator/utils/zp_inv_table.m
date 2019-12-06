function zp_inv = zp_inv_table(zp)

if(nargin < 1 || zp == 30097)
    p = mfilename('fullpath');
    p = [p(1:end-12) 'zp_linsolve_ik30097.mat'];
    zp_inv = load(p);
    zp_inv = zp_inv.ik;
else
    zp_inv = zeros(1,zp-1);
    for k = 1:zp-1
        zp_inv(k) = zp_inverse(k,zp);
    end
end
