function y = prod(x)
% MULTIPOL prod
if(size(x,1)==1); x = x'; end;
for j = 1 : size(x, 2)
    y(1,j) = multipol(1, zeros(nvars(x(1)),1));
    for i = 1 : size(x, 1)
        y(1,j) = y(1,j)*x(i,j);
    end
end