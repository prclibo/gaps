function m = monvec2matrix(mons)
% m = monvec2matrix(mons) takes a vector of monomials and returns a
% matrix by concatenating all exponents. monvec2matrix([x^2y xy]) => 
% [2 1;1 1]
m = zeros(nvars(mons(1)), length(mons));
for k = 1 : length(mons)
    if(nterms(mons(k)) ~= 1)
        error('input has to be a vector of monomials (not polynomials)');
    end
    m(:, k) = monomials(mons(k));
end