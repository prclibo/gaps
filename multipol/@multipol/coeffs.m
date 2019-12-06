function c = coeffs(p1);
% MULTIPOL/COEFFS operator
% Returns the coefficients of the polynom p1
% c = coeffs(p1);

c={p1.coeffs};
if numel(c) == 1, c = c{1}; end
