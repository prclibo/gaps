function c = rdivide(p1,p2)
% MULTIPOL/RDIVIDE Called for p1./p2 syntax. Only handles monomials.
% For polynomial division see MULTIPOL/MRDIVIDE

for i=1:numel(p2)
	if isa(p2(i),'multipol') && nterms(p2(i))~=1
		error('Only division by monomials supported');
	end
end

if isnumeric(p1)
	c = multipol(p1)./p2;
	
elseif isnumeric(p2)
	c = p1./multipol(p2);
	
elseif isa(p1,'multipol') && isa(p2,'multipol')
	[p1 p2] = eqsize(p1,p2);
	if numel(p1)==1
		c = p2;
		for i=1:numel(p2)
			c(i).coeffs = p1.coeffs./p2(i).coeffs;
			c(i).monomials = bsxfun(@minus,p1.monomials,p2(i).monomials);
		end
	elseif numel(p2)==1
		c = p1;
		for i=1:numel(p1)
			c(i).coeffs = p1(i).coeffs./p2.coeffs;
			c(i).monomials = bsxfun(@minus,p1(i).monomials,p2.monomials);
		end
	elseif all(size(p1)==size(p2))
		c = p1;
		for i=1:numel(p1)
			c(i).coeffs = p1(i).coeffs./p2(i).coeffs;
			c(i).monomials = bsxfun(@minus,p1(i).monomials,p2(i).monomials);
		end
	else
		error('Size mismatch');
	end
else
	error('Unsupported operation');
end
c = squeeze(c);