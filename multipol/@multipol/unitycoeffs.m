function pout = unityCoeffs(pin)
%pout = unityCoeffs(pin)
%returns a polynomial pout with the same monomials as pin but
%with all coefficients set to 1.
m = pin.monomials;
c = ones(1,size(m,2));
pout = multipol(c,m);
