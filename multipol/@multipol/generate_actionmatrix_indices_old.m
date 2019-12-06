function ind = generate_actionmatrix_indices_old(monoms, basisind, var)
% N.B! this function exists for compability reasons. Use
% generate_actionmatrix_indices instead.
%
% ind = generate_actionmatrix_indices_old(monoms, basisind, var) generates a
% vector of indices off-line for use in gb-solvers.

n_mons = nterms(monoms);
mons = monomials(monoms);
mm = multipol([1:n_mons], monomials(monoms));

% loop through all monomials in monoms specified by basisind
for k = 1 : length(basisind)
    mon = mons(:,basisind(k));
    
    % multiply each with var
    mon_var = multipol(1,mon)*var;
    
    % check where it ends up after multiplication
    ind(k) = coeffsof(mon_var, mm);
end