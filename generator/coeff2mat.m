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

function [C, M] = coeff2mat(coefs, mdegs, order)

if nargin<3
	order = 'grevlex';
end

if ~iscell(coefs), coefs = {coefs}; end
if ~iscell(mdegs), mdegs = {mdegs}; end

neqs = numel(coefs);
nterms = cellfun(@(x) size(x, 2), mdegs);

nt = sum(nterms);
nv = size(mdegs{1}, 1);
assert(all(cellfun(@(x) size(x, 1), mdegs) == nv));

M = zeros(nv, nt,'int32');
inds = cell(neqs, 1);
k = 0;
for i=1:neqs
	inds{i} = k + (1:nterms(i));
	M(:, k + (1:nterms(i))) = mdegs{i};
	k = k + nterms(i);
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

C = zeros(neqs, size(M,2), class(coefs{1}));
% if isa(coefs{1}, 'sym'), C = sym(C); end
for i=1:neqs
	ind = ib(inds{i});
	C(i,ind) = coefs{i};
end
