function [mod V stats] = construct_modulomatrix(C, mon, dim, P, x, method)
% [mod V stats] = construct_modulomatrix(C, mon, dim, P, x, method) constructs the matrix for
% the modulo mapping (reduction modulo I).
%
% INPUT
%   C - coefficient matrix
%   mon - vector of monomials
%   dim - expected number of solutions
%   P - index vector of possible (permissible) basis monomials
%   x - multipol action variable (set to string 'all' if not specified, this is default)
%   method - 'qr', 'std' or ('svd') (default 'qr');
%
% OUTPUT
%   mod - modulo matrix
%   V - permutation vector for 'qr', change-of-basis matrix for 'svd',
%   empty for 'std'.
%   stats - some statistics

smart_elim = false; % costs some time but increases chance of getting a
% solution. not yet working as it should.
 

if(nargin < 6), method = 'qr'; end;

implemented_methods = {'qr'};

if(~ismember(method, implemented_methods))
    error(['method: ' method ' not implemented yet']);
end

p = nvars(mon(1));
n = length(mon);
r = dim;

if(isstr(x) && strcmp(x, 'all'))
    ind = [];
    for k = 1 : p
        xk = zeros(p, 1);
        xk(k) = 1;
        xk = multipol(1, xk);
        ind_k = generate_actionmatrix_indices(mon, P, xk);
        ind = union(ind, ind_k);
    end
else
    ind = generate_actionmatrix_indices(mon, P, x);
end

% monomials to reduce
R = setdiff(ind, P);
all = 1 : n;
% excessive monomials
E = setdiff(all, union(R, P));

% reorder
ne = length(E);
nr = length(R);
np = length(P);
stats.n_to_reduce = nr;
stats.n_excessive = ne;
V = [E R P];
C = C(:, [E, R, P]);
mon = mon([E, R, P]);

% eliminate excessive monomials (done by lu in the qr paper.
% this is more general)

% the row_echelon method was used before, but the qr version is almost
% always much faster.
%[C2 k] = row_echelon(C, length(E));
if(~isempty(E))
    [qq rr ee] = qr(C(:, 1 : length(E)));
    C2 = [rr qq'*C(:, length(E) + 1 : end)];
    kk = abs(rr(1))./abs(diag(rr)) < 1e10;
    k = find(diff(kk) == -1);
    %keyboard
    if(isempty(k))
        k = length(kk);
    end
else
    C2 = C;
    k = 0;
end

% partition C into R- and P-parts
CR = C2(k + 1 : end, ne + 1 : ne + nr);
CP = C2(k + 1 : end, end - np + 1 : end);
mm = size(CR, 1);
nn = size(CR, 2) + size(CP, 2);
if(nn - mm > r)
    fprintf('rank(CR) %d size(CR) [%d %d]', rank(CR), size(CR,1), size(CR,2));
    error(['Not enough equations for that solution dimension, r: '...
        num2str(r) ', nn: ' num2str(nn) ', mm: ' num2str(mm)...
        '. nn-mm <= r does not hold. Try increasing the basis size']);
end

% eliminate R-monomials (this step is included in the lu factorization
% in the paper. qr is slightly slower but more stable).
[q2 UR2_0] = qr(CR);
CP = q2'*CP;

% select basis (qr + column pivoting)
CP2 = CP(1 : nr, :);
CP3 = CP(nr + 1 : end, :);
[q3 r3 e] = qr(CP3, 0);
CP4 = CP2(:, e(1 : end - r));
CB1 = CP2(:, e(end - r + 1 : end));
UP3 = r3(1 : np - r, 1 : np - r);
CB2 = r3(1 : np - r, end - r + 1 : end);
if(isempty(CP4)), CP4 = []; end;
if(isempty(UP3)), UP3 = []; end;
ee = [1 : ne + nr e + ne + nr];
V = V(ee);
mon = mon(ee);

%
% elimination step
%
Celim = [UR2_0(1 : nr + np - r, :) [CP4; UP3]];

% smart elim, i.e. renew check of which monomials to reduce.
% todo: write extra-smart elim by iteratively doing echelon and moving
% monomials from the permissible set to the non-permissible until convergence.
% todo: there is a bug here. fix.
if(smart_elim)
    B2 = nr + np - r + 1 : nr + np;
    ind = generate_actionmatrix_indices(mon, B2 + ne, x) - ne;
    R2 = setdiff(ind, B2);
    E2 = setdiff(1 : nr + np, union(B2, R2));
    V = V([1:ne E2+ne R2+ne B2+ne]);
    mon = mon([1:ne E2+ne R2+ne B2+ne]);
    Csub = [Celim [CB1; CB2]];
    Csub = Csub(:, [E2 R2 B2]);
    [Ce2, k] = row_echelon(Csub, length(E2));
    Celim = Ce2(k + 1 : end, length(E2) + 1 : end - r);
    CB = Ce2(:, end - r + 1 : end);
    T = - Celim \ CB(k + 1 : end, :);
    T = [zeros(k, r); T];
else
    T = - Celim \ [CB1; CB2];
% 	T = -pinv(Celim)*[CB1; CB2];
end

stats.rankdiff = size(Celim, 2) - rank(Celim);
stats.condition = cond(Celim);

% modulo matrix
mod = zeros(r, n);
mod(:, end - r + 1 : end) = eye(r);
mod(:, ne + 1 : end - r) = T';