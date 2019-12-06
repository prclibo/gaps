function p3 = subs(p1,var,p2);
% MULTIPOL/SUBS operator
% Substitution
% p3 = subs(p1,var,p2);
% p3 = p1(..,p2,..); where the monomial var is substituted for p2.
% p2 is a double or a multipol object.

if isa(p1, 'cell')
    polys = p1;
else
    polys = {p1};
end

for k = 1 : length(polys)
    p1 = polys{k};
    c = coeffsof(var, p1);
    p3{k} = p1 - c * var;
    p3{k} = p3{k} + c * p2;
end

if length(p3) == 1,
    p3 = p3{1};
end