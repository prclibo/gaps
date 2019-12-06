function r = gt(p1,p2)

[p1 p2] = eqsize(p1,p2);

m1 = p1.monomials;
m2 = p2.monomials;

% Larger if larger number of terms
if size(m1,2)>size(m2,2), r = true; return; end
if size(m1,2)<size(m2,2), r = false; return; end

% Lex(x,y,z,...) order monomials
m1 = -sortrows(-m1')';
m2 = -sortrows(-m2')';

s1 = sum(m1,1);
s2 = sum(m2,1);

for i=1:size(m1,2)
	if s1(i)>s2(i), r = true; return; end
	if s1(i)<s2(i), r = false; return; end
end

d = m1-m2;
i = find(d~=0,1);
if ~isempty(i)
	if d(i)>0, r = true;
	else r = false;
	end
	return;
end

% All monomials equal, order by coefficients
dc = p1.coeffs-p2.coeffs;
i = find(dc~=0,1);
if ~isempty(i)
	if dc(i)>0, r = true;
	else r = false;
	end
	return;
end

% p1 == p2
assert(p1==p2)
r = false;