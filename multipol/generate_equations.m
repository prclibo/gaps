function [eqs,mon] = generate_equations(eqs,degree,multmons)
% eqs = generate_equations(eqs, degree) generates an expanded set of
% equations by multiplying with all monomials up to degree.
eqs = eqsize(eqs);
if nargin<3
	p = nvars(eqs(1));
	x(1) = multipol(1, zeros(p, 1));
	for k = 1 : p
		var = zeros(p, 1);
		var(k) = 1;
		x(k + 1) = multipol(1, var);
	end
else
	x = eqsize(multmons);
end

mon = 1;
for k = 1 : degree
	mon = mon(:) * x(:)';
	mon = monvec(sum(mon(:))); % this line removes duplicate monomials.
end
eqs = eqs(:) * mon(:)';
eqs = eqs(:);
