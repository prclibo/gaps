function [sols stats] = polysolve(eqs, settings)
% [sols stats] = polysolve(eqs, settings)
%
% POLYSOLVE solves a set of polynomial equations in p variables using
% action matrix techniques. It is not guaranteed to work for all inputs.
% Sucess depends on the complexity of the equations and the parameters
% in settings.
%
% INPUT:
%   eqs             - A vector of multipol objects.
%   settings - A struct with the following
%     required fields:
%       dim - Number of solutions.
%     optional fields:
%       basis_size - Number of permissible basis mons to use. Set to 'all'
%           to use all possible permissible monomials. (default 'all').
%       method - 'qr', 'svd' or 'std' (default 'qr').
%       action_variable - Index for the action variable (i.e. =1 means
%                         use the first variable). If set to 'all'
%                         action matrices are generated for all
%                         variables and solutions are extracted from
%                         the eigenvectors using a random linear
%                         combination of action matrices. This is the
%                         safest choice. (default 'all').
%
% OUTPUT:
%   sols            - a p by dim matrix of solutions.
%   condition       - condition number of the elimination matrix. For ok accuracy
%                   should usually be < 1e8 to 1e10.

% check fields
fields = fieldnames(settings);
if(~ismember('dim', fields)), error('settings.dim required'); end
if(~ismember('basis_size', fields)), settings.basis_size = 'all'; end
if(~ismember('action_variable', fields)), settings.action_variable = 'all'; end
if(~ismember('method', fields)), settings.method = 'qr'; end

p = nvars(eqs(1));
r = settings.dim;

% define action variable
if(~strcmp(settings.action_variable, 'all'))
    p = nvars(eqs(1));
    x = zeros(p, 1);
    x(abs(settings.action_variable)) =...
        sign(settings.action_variable) * 1;
    x = multipol(1, x);
else
    x = 'all';
end

% get coeff matrix
fprintf('extracting coefficient matrix...\n');
[C mon] = polynomials2matrix(eqs);
n = length(mon);
% decide on possible basis monomials
fprintf('checking for permissible basis monomials...\n');
if(isfield(settings, 'permissible'))
    P = settings.permissible;
else
    [non_perm perm] = getPermissible(mon,...
        settings.action_variable);
    r = settings.dim; % number of solutions
    
    np = settings.basis_size;
    if(strcmp(np, 'all'))
        P = perm;
    else
        if(np > length(perm))
            error('basis_size too large');
        end
        P = perm(end - np + 1 : end);
    end
end

if(length(P) < settings.dim)
    error('not enough possible basis monomials: p < r');
end

%
% construct modulo matrix
%
fprintf('constructing modulo matrix...\n');
[mod e modstats] = construct_modulomatrix(C, mon, settings.dim,...
    P, x, settings.method);
mon = mon(e);

%
% save some statistics
%
stats.n_monomials = length(mon);
stats.n_eqs = size(C, 1);
stats.rank = rank(C);
stats.rankdiff = size(C, 2) - stats.rank;
stats.inner_rankdiff = modstats.rankdiff;
stats.condition = modstats.condition;
stats.n_permissible = length(P);
stats.n_to_reduce = modstats.n_to_reduce;
stats.n_excessive = modstats.n_excessive;
stats.reducible_range = [stats.n_excessive + 1, stats.n_monomials];
reducibles = stats.reducible_range(1) : stats.reducible_range(2);
stats.basis = mon(end - r + 1 : end);

%
% construct action matrix
%polysolve_test2.m
fprintf('constructing action matrix...\n');
m = construct_actionmatrix(mon, mod, settings.dim, P, x);

fprintf('eigen decomposing...\n');
if(isstr(x) && strcmp(x, 'all'))
    for k = 1 : p
        [vv, dk] = eig(m{p}');
        d(:, k) = diag(dk);
    end
else
    [vv, tmp] = eig(m');
end

if(iscell(m))
    mm = m;
    m = m{end};
end
if(false && strcmp(x, 'all'))
    sols = d';
else
    if 0;
        % normalize eigenvectors
        one = multipol(1, zeros(p, 1));
        onev = zeros(n, 1);
        onev(mon == one) = 1;
        disp(find(onev));
        onered = mod * onev;
        vvn = vv/diag(vv.'*onered);
        
        % compute solutions (values of all p xk on the solution set)
        sols = zeros(p, r);
        for k = 1 : p
            % define xk
            v = zeros(p, 1);
            v(k) = 1;
            xk = multipol(1, v);
            
            % find xk among the monomials and define a corresponding
            % coefficient vector vk.
            vk = zeros(n, 1);
            ind = find(mon == xk);
            vk(ind) = 1;
            if(~ismember(ind, reducibles))
                % this should never happen as the algorithm is currently
                % written. (unless perhaps in the unlikely case of more
                % variables than solutions).
                error(['action matrix computation ok, but could not compute value of variable ', k]);
            end
            
            % reduce to the quotient space if necessary
            xred = mod * vk;
            sols(k, :) = xred' * vvn; % xk is a linear combination of the solving basis elements.
        end
        
        
        
        
    else
        % normalize eigenvectors, version 2
        version2 = 1;
        version1 = 0;
        %keyboard
        if version1,
            MM = zeros(p,r);
            for kkk = 1:r;
                MM(:,kkk) = monomials(mon( n-r+kkk));
            end
            MMe = [MM;ones(1,r)];
        end;
        if version2,
            okmon = find(sum(mod,1)~=0);
            MM = zeros(p,length(okmon));
            for kkk = 1:length(okmon);
                MM(:,kkk) = monomials(mon(okmon(kkk)));
            end
            MMe = [MM;ones(1,length(okmon))];
        end
        %z = null(MM);
        %z = z./(repmat(sum(z,1),r,1));
        %blubb =  abs(exp(-z'*log(abs(vv))));
        %vvn = vv*diag(blubb(1,:));
        %blubb2 =  abs(exp(-z'*log(abs(vvn))));
        
        % ber???kna l???sningarna
        if 0,
            asol = (MMe')\log(abs(vv));
            err1 = log(abs(vv)) - (MMe')*asol;
        end;
        % H???r kan man f??? en indikation p??? om det g???tt bra eller inte.
        
        % En annan variant ???r att leta efter submatriser i MMe
        % som har determinant 1.
        bestN = Inf;
        bestbb = NaN;
        for kkk = 1:100;
            bb = randperm(size(MMe,2));
            thisN = round(abs(det(MMe(:,bb(1:(p+1))))));
            if thisN>0.1,
                if thisN<bestN,
                    bestN = thisN
                    bestbb = bb(1:(p+1));
                end
            end
        end;
        
        if version1,
            if bestN >= 1,
                xsol = exp((MMe(:,bestbb)')\log(vv(bestbb,:)));
                err = real((MMe')*log(xsol) - log(vv ));
                ok = log10(sqrt(sum(abs(err))));
                oki = ok<-1;
                stats.ok = ok;
                err2 = (MMe')*log(xsol) - log(vv );
                fixx = round(imag(err2)/(2*pi));
                xsol2 = exp( (MMe')\(log(vv)+fixx*2*pi*i) );
            end
        end;
        if version2,
            if bestN >= 1,
                vv = mod(:,okmon)'*vv;
                xsol = exp((MMe(:,bestbb)')\log(vv(bestbb,:)));
                err = real((MMe')*log(xsol) - log(vv ));
                ok = log10(sqrt(sum(abs(err))));
                oki = ok<-1;
                stats.ok = ok;
                err2 = (MMe')*log(xsol) - log(vv );
                fixx = round(imag(err2)/(2*pi));
                xsol2 = exp( (MMe')\(log(vv)+fixx*2*pi*i) );
            end
        end;
        
        % N???sta m???l ???r att hitta en heltalsmatris fixx s??? att
        % avv + fixx ligger i kolumnrummet f???r MMe'
        %disp('Nu ska vi se hur det gick');
        %keyboard;
        sols = xsol2(1:(end-1),:);
        
    end;
    
end
p=0;

