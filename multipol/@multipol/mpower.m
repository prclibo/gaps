function p2 = mpower(p1,a)
% MULTIPOL/MPOWER operator
% Power of a multipol.
% p2 = p1^a;

if ~isnumeric(a) || any(round(a)~=a)
	error('Exponent must be an integer');
end

if numel(a)==1
	p2 = p1;
	for i=1:numel(p1)

		p2(i) = multipol(1,zeros(size(p1(i).monomials,1),1), 'vars', p1.vars);
		for ii = 1:a
			p2(i) = p2(i)*p1(i);
		end

	end
	
elseif numel(p1)==1
	
	for i=numel(a):-1:1
		p2(i) = p1^a(i);
	end
	p2 = reshape(p2,size(a));
	
elseif all(size(p1)==size(a))
	
	p2 = p1;
	for i=1:numel(p1)
		p2(i) = p1(i)^a(i);
	end
	
else
	error('Dimension mismatch');
end