function coeffs = coeffsof(monoms, polynomial)
%coeffs = coeffsOf(monoms, polynomial)
%returns a vector of coefficients of the monomials in "monoms" taken from
%"polynomial". Both monoms and polynomial are multipol objects.

nMonoms = nterms(monoms);
nTerms = nterms(polynomial);

mp = polynomial.monomials;
m = monoms.monomials;
m = repmat(m,[1,1,nTerms]);
m = permute(m, [1 3 2]);
mp = repmat(mp,[1,1,nMonoms]);

comp = (m==mp)*1;
comp = prod(comp, 1);
if(length(size(comp))==3)
    % we are extracting more than one coefficient
    ccomp(:,:) = comp(1,:,:);
    
    missingMons = (sum(ccomp, 1) == 0);
    ccomp = [missingMons; ccomp];
    [I J] = find(ccomp);
else
    % we are extracting only one coefficient
    if(sum(comp)==0)
        I = 1;
    else
        I = find(comp) + 1;
    end
end

c = polynomial.coeffs;
c = [0 c];

coeffs = c(I);
