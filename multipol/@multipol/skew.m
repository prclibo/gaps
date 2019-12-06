function A=skew(a);
% A=skew(a) - returns skew matrix y such that a x v = Av for any v
%
m = monomials(a(1));
n = size(m,1);
zero = multipol(0,zeros(n,1));
A=[zero,-a(3),a(2);a(3),zero,-a(1);-a(2),a(1),zero];
