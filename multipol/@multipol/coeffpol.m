function coeffpolynomial = coeffpol(polynomial, variable)
% coeffpolynomial = collect(polynomial, variable)
% well, this function is a bit weird...

if(size(variable)>1); error('coeffpol: variable has to be a monomial');end;
mon = monomials(variable);
mons = monomials(polynomial);
cfs = coeffs(polynomial);
if(length(find(mon)) > 1); error('coeffpol: variable has to be a single variable'); end;

ind = find(mon);
degree = mon(ind);

ind2 = find(mons(ind,:)==degree);

retmons = mons(:,ind2);
retmons(ind,:) = retmons(ind,:) - degree;

retcoeffs = cfs(ind2);

coeffpolynomial = multipol(retcoeffs, retmons);