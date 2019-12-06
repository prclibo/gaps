function [A,jb,row_ind] = zp_rref(A,zp)
zp_inv = zp_inv_table(zp);
A = round(A);
A = mod(A,zp);


% Adapted from MATLABs rref
[m,n] = size(A);
i = 1;
j = 1;
jb = [];
row_ind = 1:size(A,1);
while (i <= m) && (j <= n)
    % Find value and index of largest element in the remainder of column j.
    
    k = find(A(i:m,j),1);    
    if isempty(k)
        j = j + 1;
        continue;
    end    
    k = k + i - 1;
    
    % Remember column index
    jb = [jb j];

    % Swap i-th and k-th rows.
    A([i k],j:n) = A([k i],j:n);
    row_ind([i k]) = row_ind([k i]);

    % Divide the pivot row by the pivot element.
    A(i,j:n) = mod(A(i,j:n) * zp_inv(A(i,j)), zp);

    % Subtract multiples of the pivot row from all the other rows.
    for k = [1:i-1 i+1:m]
        A(k,j:n) = mod(A(k,j:n) - mod(A(k,j)*A(i,j:n),zp),zp);
    end
    i = i + 1;
    j = j + 1;
end