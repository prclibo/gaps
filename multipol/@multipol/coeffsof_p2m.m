function coeffs = coeffsof_p2m(monoms, polynomial)
%coeffs = coeffsOf(monoms, polynomial)
%returns a vector of coefficients of the monomials in "monoms" taken from
%"polynomial". Both monoms and polynomial are multipol objects.

nMonoms = nterms(monoms);
nTerms = nterms(polynomial);

monom_ind = zeros(1,nTerms);

for ii = 1:nTerms;
    monom_ind(ii) = find(sum(...
        abs(monoms.monomials-polynomial.monomials(:,ii)*...
        ones(1,nMonoms)), 1)==0);
end

coeffs = zeros(1,nMonoms);
coeffs(monom_ind) = polynomial.coeffs;
