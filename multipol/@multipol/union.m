function p = union(p1,p2)
% MULTIPOL/UNION Returns the union of the polynomials in sets p1 and p2
p = unique([p1(:); p2(:)]);
