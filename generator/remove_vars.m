function [ eqs ] = remove_vars( eqs, vars )

sz = size(eqs);
eqs = eqs(:);
[cc,mm] = polynomials2matrix(eqs(:));
mmv = monvec2matrix(mm);
if any(any(mmv(vars,:)))
    % We drop any monomials which contain the variables.
    ind = any(mmv(vars,:),1);
    cc(:,ind) = [];
    mmv(:,ind) = [];
    warning('Removing nonzero exponents.')
end
mmv(vars,:) = [];

%eqs = cc*mm;
eqs = m2p(cc,mmv);
eqs = reshape(eqs,sz(1),sz(2));

end

