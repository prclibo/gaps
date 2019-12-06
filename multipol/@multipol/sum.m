function S=sum(X,dim)
% SUM Sum of elements.
%    S = SUM(X) is the sum of the elements of the vector X. If
%    X is a matrix, S is a row vector with the sum over each
%    column. 
% 
%    S = SUM(X,DIM) sums along the dimension DIM. 

% If dimension is given transpose X if neccessery and then cal sum again.
if nargin == 2
    if dim == 2
        S = sum(X')';
    else
        S = sum(X);
    end
else
    % If it is an vector
    if  min(size(X))==1
        S = multipol();
        for ii = 1:length(X);
            S = S+X(ii);
        end
    else
        % If X is a matrix sum over the columns
        for ii = 1:size(X,2);
            S(ii) = sum(X(:,ii));
        end
    end
end
