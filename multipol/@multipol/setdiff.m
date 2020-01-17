function [p ia] = setdiff(eq1,eq2)
% MULTIPOL/SETDIFF Returns the polynomials in eq1 not also present in eq2
% ia are indices such that p = eq1(ia)


% Naive implementation, should use sorting instead.
% 
% ia = true(size(eq1));
% for i=1:numel(eq1)
% 	for j=1:numel(eq2)
% 		if eq1(i)==eq2(j)
% 			ia(i) = false;
% 			break;
% 		end
% 	end
% end
% 
% p = eq1(ia);
% ia = find(ia);


% With sorting
%
% [eq1 ind] = sortpolys(eq1);
% eq2 = sortpolys(eq2);
% 
% ia = true(size(eq1));
% j = 1;
% for i=1:numel(eq1)
% 	while j<=numel(eq2) && eq1(i)<eq2(j), j=j+1; end
% 	if j>numel(eq2), break; end
% 	if eq1(i)==eq2(j), ia(i) = false; end
% end
% 
% p = eq1(ia);
% ia = ind(ia);

% assert(all(ia==ia1))

% % Hashing and using builtin setdiff, much faster
% [eq1 eq2] = eqsize(eq1,eq2);
% 
% c1 = cell(size(eq1));
% c2 = cell(size(eq2));
% for i=1:numel(eq1)
% 	c1{i} = hash(eq1(i));
% end
% for i=1:numel(eq2)
% 	c2{i} = hash(eq2(i));
% end
% 
% [~,ia] = setdiff(c1,c2,'stable');
% ia = ia(:)';
% p = eq1(ia);

[eq1 eq2] = eqsize(eq1,eq2);
[~,ia] = setdiff(string(eq1), string(eq2),'stable');
ia = ia(:)';
p = eq1(ia);


