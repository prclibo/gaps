%POLYNOMIALS2MATRIX Decomposes a system of multipol equations into a
%	coefficient matrix and a vector of monomials.
%
%	[C MON] = POLYNOMIALS2MATRIX(P,[order]) returns a vector MON of the
%	unique monomials in the polynomials P and a matrix C so that P = C*MON.
%	Possible monomial orders are
%		'same'   : don't change the order
%		'plex'   : sort in pure lexicographical order
%		'grlex'  : sort in graded lex order
%		'grevlex': sort in graded reverse lex order (default)

function [C mon] = polynomials2matrix(p,order)

if nargin<2
	order = 'grevlex';
end

p = eqsize(p);
nt = sum(nterms(p)); % This is the slow part...
nv = nvars(p(1));
unk_vars = vars(p(1));

M = zeros(nv,nt,'int32');
inds = cell(numel(p),1);
k = 1;
for i=1:numel(p)
	inds{i} = k:k+nterms(p(i))-1;
	M(:,k:k+nterms(p(i))-1) = monomials(p(i));
	k = k + nterms(p(i));
end

switch order
	case 'same'
		[~,ia,ib] = unique(M','rows','stable');
	case 'plex'
		[~,ia,ib] = unique(-M','rows');
	case 'grlex'
		[~,ia,ib] = unique([-sum(M,1)' -M'],'rows');
	case 'grevlex'
		[~,ia,ib] = unique([-sum(M,1)' fliplr(M')],'rows');
	otherwise
		error('Unknown monomial order');
end
M = double(M(:,ia));

for i=size(M,2):-1:1
	mon(i,1) = multipol(1, M(:,i), 'vars', unk_vars);
end

%C = zeros(numel(p),size(M,1));
C = zeros(numel(p),size(M,2)); % bugfix Viktor 20160816
if isa(coeffs(p(1)), 'sym'), C = sym(C); end
for i=1:numel(p)
	ind = ib(inds{i});
	C(i,ind) = coeffs(p(i));
end
