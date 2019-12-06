function e = eq(p1,p2)

if(numel(p1) > 1 && numel(p2) == 1)
	e = false(size(p1));
	for i=1:numel(p1)
		e(i) = p1(i) == p2;
	end
	return;
end
if(numel(p1) == 1 && numel(p2) > 1)
	e = false(size(p2));
	for i=1:numel(p2)
		e(i) = p2(i) == p1;
	end
	return;
end
if(numel(p1) > 1 && numel(p2) > 1)
	if(~(all(size(p1) == size(p2)))), error('Inputs must be of equal length'); end
	e = false(size(p1));
	for i = 1 : size(p1, 1)
		for j = 1 : size(p1, 2)
			e(i,j) = p1(i,j) == p2(i,j);
		end
	end
	return;
end

if isnumeric(p1)
	p1 = multipol(p1);
end
if isnumeric(p2)
	p2 = multipol(p2);
end
[p1 p2] = eqsize(p1,p2);
e = false;
if isnumeric(p1)
	if any(bsxfun(@times,p2.coeffs,p2.monomials)~=0), return; end
	if p2.coeffs(sum(p2.monomials,1)==0)~=p1, return; end % there should be at most one constant monomial, if not something will break.
	e = true;
	return;
end
if isnumeric(p2)
	if any(bsxfun(@times,p1.coeffs,p1.monomials)~=0), return; end
	if p1.coeffs(sum(p1.monomials,1)==0)~=p2, return; end
	e = true;
	return;
end

c1 = p1.coeffs; c2 = p2.coeffs;
if numel(c1) ~= numel(c2), return; end
if isempty(c1) && isempty(c2), e = true; return; end
if ~all(c1==c2), return; end

m1 = p1.monomials; m2 = p2.monomials;
if size(m1,1) ~= size(m2,1), return; end
if ~all(m1==m2), return; end
% if(sum(abs(coeffs(p1) - coeffs(p2)))~=0), return; end
% if(sum(sum(abs(monomials(p1) - monomials(p2)))) ~= 0), return; end
e = true;