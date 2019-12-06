function ind = generate_actionmatrix_indices(monoms, basisind, var)
% ind = generate_actionmatrix_indices(monoms, basisind, var) generates a
% vector of indices off-line for use in gb-solvers.

if(length(monoms) == 1)
    ind = generate_actionmatrix_indices_old(monoms, basisind, var);
    return;
end

n_mons = length(monoms);
mons = monvec2matrix(monoms);
mm = multipol(1:n_mons, mons);

% loop through all monomials in monoms specified by basisind
for k = 1 : length(basisind)
    mon = mons(:,basisind(k));
    
    % multiply each with var
    mon_var = multipol(1,mon)*var;
    
    % check where it ends up after multiplication
    ind(k) = coeffsof(mon_var, mm);
end