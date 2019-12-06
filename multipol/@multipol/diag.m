function a=diag(A)
% DIAG Diagonal matrices and diagonals of a matrix.
%    DIAG(V) when V is a vector with N components is a square matrix
%    of order N.
% 
%    DIAG(X) when X is a matrix is a column vector formed from
%    the elements of the diagonal of X.
% 

% If input is a vector
if min(size(A))==1
    for ii = 1:length(A);
        a(ii) = multipol();
    end
    a = a'*a;
    for ii = 1:length(a);
        a(ii,ii) = A(ii);
    end
else
    for ii = 1:min(size(A))
        a(ii) = A(ii,ii);
    end
    a=a';
end

