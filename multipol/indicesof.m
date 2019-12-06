function ind = indicesof(x, y)
% ind = indicesof(x, y) finds the positions of the numbers y in the vector
% x. It is assumed that the vector x contains only unique elements.

if(length(unique(x)) ~= length(x))
    error('x must contain only unique elements');
end

if(~isempty(setdiff(y, x)))
    error('some numbers in y cannot be found in x');
end

% repair holes
extra = setdiff(1 : max(x), x);
x = [x extra];

% sort
[x_sorted ind] = sort(x);

% find what we are looking for
ind = ind(y);