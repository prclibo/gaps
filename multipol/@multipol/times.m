function p3 = times(p1,p2);
% MULTIPOL/TIMES operator

if size(p1) ~= size(p2)
    error('Matrix dimensions must agree.');
end

for i = 1 : size(p1, 1)
    for j = 1 : size(p1, 2)
        p3(i, j) = p1(i,j) * p2(i,j);
    end
end