function p3 = mtimes(p1,p2)
% MULTIPOL/MTIMES operator
% Multiplication of two polynomials
% p3 = mtimes(p1,p2);
% p3 = p1*p2;

if isempty(p1) || isempty(p2)
	p3 = [];
	return;
end

if numel(p1)==1 && numel(p2)==1
	% Multiplication with scalar
	if isnumeric(p1)
		p3 = multipol(p2.coeffs*p1,p2.monomials, 'vars', p2.vars);
		p3 = squeeze(p3);
		
	elseif isnumeric(p2)
		p3 = multipol(p1.coeffs*p2,p1.monomials, 'vars', p1.vars);
		p3 = squeeze(p3);
		
	else % Multiplication between two multipols
		if isempty(p1.monomials)
            if size(p1.monomials,1)~=size(p2.monomials,1)
                [p1, p2] = eqsize(p1,p2);
            end
			p3 = multipol(p1.coeffs*p2.coeffs, 'vars', p2.vars);
		else
			n1 = numel(p1.coeffs);
			n2 = numel(p2.coeffs);
			coeffs = kron(p1.coeffs,p2.coeffs);
			monomials = p1.monomials(:,ceil((1:n1*n2)/n2)) + ...
                        p2.monomials(:,repmat(1:n2,1,n1));
			p3 = multipol(coeffs,monomials, 'vars', p1.vars);
		end
		p3 = sort(p3);
		p3 = squeeze(p3);
	end
	
elseif issparse(p1) % Multiplication with sparse coefficient matrix
	
	[m1, n1] = size(p1);
	[m2, n2] = size(p2);
	if n1~=m2
		error('Inner matrix dimensions must agree');
	end
	p3(m1,n2) = multipol;
	for r=1:m1
		[~,cols,v] = find(p1(r,:));
		for c=1:n2
			p3(r,c) = v*p2(cols,c);
		end
	end
	
else % Multiplication between two multipol matrices
	
	if numel(p1)==1
		for i=numel(p2):-1:1
			p3(i) = p1*p2(i);
		end
		p3 = reshape(p3,size(p2));
		
	elseif numel(p2)==1
		for i=numel(p1):-1:1
			p3(i) = p1(i)*p2;
		end
		p3 = reshape(p3,size(p1));
		
	else
		[m1 n1] = size(p1);
		[m2 n2] = size(p2);
		
		if n1~=m2
			error('Inner matrix dimensions must agree');
		end
		
		p3(m1,n2) = multipol;
		for row=1:m1
			for col=1:n2
				for i=1:n1
					p3(row,col) = p3(row,col)+p1(row,i)*p2(i,col);
				end
			end
		end
	end
	
end
