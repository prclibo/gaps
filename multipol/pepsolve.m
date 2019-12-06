function [sol info] = pepsolve(eq,evar,num_eqs,varargin)

verb = 1;
method = 'none';
order = 'grevlex';
redmat = 1;
visualize = 0;

for i=1:2:length(varargin)
	switch varargin{i}
		case 'verb', verb = varargin{i+1};
		case 'method', method = varargin{i+1};
		case 'order', order = varargin{i+1};
		case 'R', redmat = varargin{i+1};
		case 'vis', visualize = varargin{i+1};
		otherwise, error('Parameter parse error');
	end
end



eq = eq(:);
[C M] = polynomials2matrix(eq,order);

fprintf('%u original equations in %u monomials\n',size(C));

% Row reduction of the system may or may not be beneficial, turned off by default
% The monomial 'order' parameter has no effect in this case
switch method
	case 'qr'
		[~,r,p] = qr(C,'vector');
		k = nnz(abs(diag(r)./max(abs(diag(r))))>1e-10); % ~rank(r)
		C2 = r(1:k,:);
		M = M(p);
	case 'qrs'
		r = qr(sparse(C),0);
		C2 = full(r);
	case 'lu'
		[~,u,~,q] = lu(sparse(C),'vector');
% 		ind = abs(diag(u)./max(abs(diag(u))))>1e-10;
		C2 = full(u(:,:));
		M = M(q);
	case 'rref'
		C2 = rref(C);
		k = nnz(sum(abs(C2),2)>1e-12);
		C2 = C2(1:k,:);
	case 'lu_rref'
		[~,u,~,q] = lu(sparse(C),'vector');
		M = M(q);
		k = nnz(abs(diag(u)./max(abs(diag(u))))>1e-10);
% 		[~,ind] = sort(abs(diag(u)),'descend');
% 		u(:,
		C2 = full([speye(k) u(:,1:k)\u(:,k+1:end)]);
		
	case 'none'
		C2 = C;
% 		k = size(C,1);
	otherwise
		error('Unknown method');
end

keep = max(abs(C2),[],2)>1e-14; % Remove null equations
C2 = C2(keep,:);

if strcmp(num_eqs,'all'), num_eqs = size(C2,1); end
if num_eqs>size(C2,1), error('Too many equations requested'); end
ind = size(C2,1)-num_eqs+1:size(C2,1);

% ind = randperm(size(C2,1),num_eqs);
% ind = 1:num_eqs;
if visualize
	imagesc(C2(ind,:)~=0);
end

[A B] = eqs2mat(C2(ind,:),M,evar,verb);

A2 = cell(size(A));


if isnumeric(redmat)
	R = A{redmat};
elseif strcmpi(redmat,'rand')
	R = rand(size(A{1}));
elseif strcmpi(redmat,'eye')
	R = eye(size(A{1}));
else
	error('Unknown reduction matrix option');
end
for i=1:length(A)
	A2{i} = R'*A{i};
end
disp(cellfun(@rank,A2)')

tic
[X E] = polyeig(A2{:});
toc

[~,ind] = sortrows([real(E) imag(E)]);
evdiff = abs(diff(E(ind)));
if any(evdiff<1e-6)
	warning('There appear to be eigenvalues of multiplicity > 1, min. diff. %g',min(evdiff));
end





n = nvars(eq(1));
const = B==1;
if nnz(const)==0
	warning('Constant term not in basis');
	good = ~(isinf(E) | isnan(E));
else
	good = ~(isinf(E) | isnan(E) | abs(X(const,:))'<1e-12);
	X = bsxfun(@rdivide,X,X(const,:));
end
E = E(good); X = X(:,good);

% Find variable monomials and assemble solution
[~,bm] = polynomials2matrix(B,'same');
bm = monvec2matrix(bm);
sol_basis = [];
sol = [];
for i=1:n
	ind = find(sum(bm,1)==bm(i,:) & sum(bm,1)>0); % Find smallest power of x_i
	[~,j] = min(bm(i,ind));
	ind = ind(j);
	if isempty(ind)
		if i==evar
			sol = [sol; E.']; %#ok
			sol_basis = [sol_basis; multipol(1,[zeros(evar-1,1); 1; zeros(n-evar,1)])]; %#ok
		else
			warning('No power of variable %u in basis',i);
		end
	else
		sol = [sol; X(ind,:)]; %#ok
		sol_basis = [sol_basis; multipol(1,bm(:,ind))]; %#ok
	end
end
info = struct('sol_basis',sol_basis,'B',B,'ev',E,'X',X);
info.A2 = A2;
info.A = A;


% Am = A2{1};
% z = multipol(1,1);
% for i=1:length(A)-1
% 	Am = Am + A2{i+1}*z^i;
% end
% info.Am = Am;

if size(sol,1)==n
% 	fprintf('Found %u solutions\n',nnz(mean(abs(evaluate(eq,sol)))<1e-6));
end

end


function [A Mred] = eqs2mat(C,M,evar,verb)

neq = size(C,1);
mind = sum(abs(C),1)>0;
M = M(mind);
C = C(:,mind);

if verb, fprintf('%u equations, %u equation monomials after reduction\n',neq,length(M)); end
if verb>1, disp(sym(M).'); end

mm = monvec2matrix(M);
evar_deg = max(mm(evar,:));
mmr = mm;
mmr(evar,:) = 0;
[mmr,~,ic] = unique(mmr','rows');
mmr = mmr';


if verb, fprintf('%u remaining monomials after hiding variable %u\n',size(mmr,2),evar); end
for i=size(mmr,2):-1:1
	Mred(i,1) = multipol(1,mmr(:,i));
end
if verb>1
	disp(sym(Mred).');
end
if verb, fprintf('Matrix polynomial degree: %u\n',evar_deg); end

m = size(mmr,2); % number of monomials in basis

A = cell(evar_deg+1,1);
for deg=0:evar_deg
	A{deg+1} = zeros(neq,m);
	for j=1:neq%min(m,neq)
		ind = mm(evar,:)==deg & C(j,:)~=0; % terms in equation j of degree i in the elimination variable
		A{deg+1}(j,ic(ind)) = C(j,ind);
	end
end
if verb
	disp('Matrix ranks:')
	disp(cellfun(@rank,A)')
end

end

