function p = sortmons(p,order)
% MULTIPOL/SORTMONS operator
% Assumes p is an array of monomials and sorts them.

if nargin<2
	order = 'grevlex';
end

[~,p] = polynomials2matrix(sum(p),order);