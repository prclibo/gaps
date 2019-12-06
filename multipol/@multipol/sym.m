function s = sym(mp)
% MULTIPOL/SYM
% Convert polynomial expression to Matlab symbolic toolbox object
% if nargin<2
% 	xyzw = true;
% end

% s = sym(zeros(size(mp)));
% n = 0;
% for k=1:numel(mp)
% 	n = max(n,nvars(mp(k)));
% end
% % If few variables and xyzw==true, use x, y, z, w instead of numbered x's.
% names = {'x','y','z','w'};
% 
% for k=1:numel(mp)
% 	
% 	if isempty(mp(k).monomials)
% 		s(k) = mp(k).coeffs;
% 	else
% 		clear x;
% 		for i=nvars(mp(k)):-1:1
% 			if n<=4 && xyzw, name = names{i}; else name = sprintf('x%u',i); end
% 			x(i) = sym(name);
% 		end
% 		c = coeffs(mp(k));
% 		a = monomials(mp(k))';
% 		for i=1:length(c)
% 			s(k) = s(k) + c(i)*prod(x.^a(i,:));
% 		end
% 	end
% 	
% end

% Use CHAR conversion routine first, should be faster

% s = sym(char(mp,xyzw));
s = sym(zeros(size(mp)));
for i = 1:numel(mp)
    monos = prod(mp(i).vars(:).^mp(i).monomials);
    s(i) = mp(i).coeffs * monos(:);
end