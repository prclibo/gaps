function ret = zp_eval(mp, unk_zp, p)

if isstruct(unk_zp)
    in_unk = [unk_zp.var];
    mp_unk = vars(mp);
    [~, locb] = ismember(in_unk, mp_unk);
    assert(all(locb > 0));
    unk_zp = [unk_zp.subs];
    unk_zp = mod(int64(unk_zp(locb)), p);
end

coeff_zp = int64(coeffs(mp));
unk_zp = int64(unk_zp(:));
mdegs = int32(monomials(mp));
assert(size(mdegs, 1) == numel(unk_zp));

power_zp = powermod(unk_zp, mdegs, p);
mono_zp = ones(1, numel(coeff_zp), 'int64');
for i = 1:numel(unk_zp)
    mono_zp = mod(mono_zp .* power_zp(i, :), p);
end


ret = mod(coeff_zp(:) .* mono_zp(:), p);
ret = mod(sum(ret), p);

end