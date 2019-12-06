function coefs = collect_coeffs(mp, mdegs)

[~, ia, ib] = intersect(mdegs.', monomials(mp).', 'rows');

all_coefs = coeffs(mp);
coefs = zeros(size(all_coefs), class(all_coefs));
% if isa(coeffs(mp), 'sym'), coefs = sym(coefs); end

coefs(ia) = all_coefs(ib);