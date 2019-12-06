function p = squeezevars(p)
% p = squeezevars(p) removes any variables from p which are not used.
iscell = true;
if ~isa(p, 'cell'),
    iscell = false;
    pp = p;
    p = cell(length(p), 1);
    
    for k = 1 : length(p)
        p{k} = pp(k);
    end
end

for k = 1 : length(p)
    m = monomials(p{k});
    c = coeffs(p{k});
    remove_ind = find(sum(m, 2) == 0);
    m(remove_ind, :) = [];
    p{k} = multipol(c, m);
end

if(~iscell)
    pp = p;
    clear p;
    for k = 1 : length(pp)
        p(k) = pp{k};
    end
elseif length(p) == 1,
    p = p{1};
end