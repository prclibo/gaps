function [p ia ib] = intersect(eq1,eq2)
% MULTIPOL/INTERSECT Returns the intersection of the sets of polynomials.

% Naive implementation
% ia = false(1,numel(eq1));
% ib = false(1,numel(eq2));
% for i=1:numel(eq1)
% 	for j=find(~ib)
% 		if eq1(i)==eq2(j)
% 			ia(i) = true;
% 			ib(j) = true;
% 			break;
% 		end
% 	end
% end
% 
% p1 = eq1(ia);
% ia1 = find(ia);
% ib1 = find(ib);

% With sorting
% [eq1 ind1] = sortpolys(eq1);
% [eq2 ind2] = sortpolys(eq2);
% 
% ia = false(1,numel(eq1));
% ib = false(1,numel(eq2));
% j = 1;
% for i=1:numel(eq1)
% 	while j<=numel(eq2) && eq1(i)<eq2(j), j=j+1; end
% 	if j>numel(eq2), break; end
% 	if eq1(i)==eq2(j), ia(i) = true; ib(j) = true; end
% end
% 
% p = eq1(ia);
% ia = ind1(ia);
% ib = ind2(ib);

% assert(all(sortpolys(p1)==p))
% assert(all(ia==ia1))
% assert(all(ib==ib1))

% Hashing and using builtin, this is much faster
[eq1 eq2] = eqsize(eq1,eq2);

c1 = cell(size(eq1));
c2 = cell(size(eq2));
for i=1:numel(eq1)
	c1{i} = hash(eq1(i));
end
for i=1:numel(eq2)
	c2{i} = hash(eq2(i));
end

[~,ia,ib] = intersect(c1,c2,'stable');
ia = ia(:)';
ib = ib(:)';
p = eq1(ia);