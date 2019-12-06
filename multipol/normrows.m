function A = normrows(A)
% A = normrows(A) normalizes the rows of A to length 1.
D = diag(1./sqrt(sum(A'.^2)));
A = D*A;
