function p3 = plus(p1,p2)
% MULTIPOL/PLUS operator
% Addition of two polynomials
% p3 = plus(p1,p2);
% p3 = p1+p2;
if ~all(size(p1)==size(p2)) && ~(numel(p1)==1 || numel(p2)==1)
	error('Arrays must be same size or one must be a scalar');
end

if numel(p1)==1 && numel(p2)>1
	p3 = p2;
	for i=1:numel(p2)
		p3(i) = p1+p2(i);
	end
	return;
end

if numel(p2)==1 && numel(p1)>1
	p3 = p1;
	for i=1:numel(p1)
		p3(i) = p2+p1(i);
	end
	return;
end

for i=numel(p1):-1:1
	if isnumeric(p1(i))
		% Check for constant term
		ind = find(sum(abs(p2(i).monomials),1)==0);
		if ~isempty(p2(i).monomials) && ~isempty(ind)
			p3(i) = p2(i);
			p3(i).coeffs(ind) = p3(i).coeffs(ind)+p1(i);
		else
			p3(i) = multipol([p1(i) p2(i).coeffs],...
                [zeros(size(p2(i).monomials,1),1) p2(i).monomials],...
                'vars', vars(p2(i)));
		end
		
	elseif isnumeric(p2(i))
		% Check for constant term
		ind = find(sum(abs(p1(i).monomials),1)==0);
		if ~isempty(p1(i).monomials) && ~isempty(ind)
			p3(i) = p1(i);
			p3(i).coeffs(ind) = p3(i).coeffs(ind)+p2(i);
		else
			p3(i) = multipol([p2(i) p1(i).coeffs],...
                [zeros(size(p1(i).monomials,1),1) p1(i).monomials],...
                'vars', vars(p1(i)));
		end
		
	else % Should be two multipols
		if size(p1(i).monomials,1)~=size(p2(i).monomials,1)
			[p1(i) p2(i)] = eqsize(p1(i),p2(i));
		end
		p3(i) = p1(i);
		p3(i).monomials = [p1(i).monomials p2(i).monomials];
		p3(i).coeffs = [p1(i).coeffs p2(i).coeffs];
	end
	p3(i) = squeeze(sort(p3(i)));
end
p3 = reshape(p3,size(p1));