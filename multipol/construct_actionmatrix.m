function m = construct_actionmatrix(mon, mod, dim, P, x, ind)
% m = construct_actionmatrix(C, mon, dim, P, x) constructs the action matrix
% for multiplication by x.
%
% INPUT
%   mod - modulo matrix
%   mon - vector of monomials
%   dim - size of the action matrix (dim x dim)
%   P - index vector of permissible basis monomials
%   x - action variable (set x to 'all' to generate a random linear
%           combination of action matrices for all variables)
%   ind - optional pre-calculated index vector.
% OUTPUT
%   m - action matrix
%

if nargin<6,
    hasind = 0;
else
    hasind = 1;
end;

p = nvars(mon(1));
n = length(mon);
r = dim; 

if(isstr(x) && strcmp(x, 'all'))
    mm = zeros(r, r, p);
    m = {};
    for k = 1 : p
        xk = zeros(p, 1);
        xk(k) = 1;
        xk = multipol(1, xk);
        mm(:,:,k) = construct_actionmatrix(mon, mod, dim, P, xk);
        m{end + 1} = mm(:,:,k);
    end
    mm = reshape(mm, [r*r, k]);
    mm = mm * randn(k, 1);
    mm = reshape(mm, [r, r]);
    m{end+1} = mm;
    return;
end

% find action indices
if ~hasind,
  ind = generate_actionmatrix_indices(mon, n - r + 1 : n, x);
  %keyboard;
end;

% construct m
M = zeros(n, r);
M(ind, :) = eye(r);
m = mod * M;

% remove NaNs to prevent crash
if(any(isnan(m))) 
    warning('NaNs encountered during construction of the action matrix');
    m(isnan(m)) = 0;
end
