function str = cg_setup_action_matrix(solv,template,opt)
str = [];

% str = [str sprintf('C0 = C(:,1:%d);\n', template.C_sz(2)-length(template.basis))];
% str = [str sprintf('C1 = C(:,%d:end);\n', template.C_sz(2)-length(template.basis)+1)];

if opt.generalized_eigenvalue_solver
    str = [str generalized_eigenvalue_solver(solv,template,opt)];
    return;
end

% check if the template is square
if template.C_sz(1) == template.C_sz(2)-length(template.basis)
    if strcmp(opt.linear_solver,'backslash')
        str = [str sprintf('C1 = C0 \\ C1;\n')];
    elseif strcmp(opt.linear_solver,'qr')
        error('NYI: qr for square template');
    elseif strcmp(opt.linear_solver,'pinv')
        str = [str sprintf('C1 = pinv(C0)*C1;\n')];
    end
    str = [str sprintf('RR = [-C1(end-%d:end,:);eye(%d)];\n',length(template.reducible)-1,length(template.basis))];
else
    % solve in least squares sense (Non-square template)
    b = zeros(template.C_sz(2)-length(template.basis),length(template.reducible));
    b(end-length(template.reducible)+1:end,:) = -eye(length(template.reducible));
    b_ind = find(b);
    b_sz = size(b);
    
    str = [str cg_format_vector('b_ind',b_ind,opt.max_str_len)];
    if opt.sparse_template
        str = [str sprintf('b = sparse(%d,%d);\n',b_sz(1),b_sz(2))];
    else
        str = [str sprintf('b = zeros(%d,%d);\n',b_sz(1),b_sz(2))];
    end
    str = [str sprintf('b(b_ind) = -1;\n')];

    if strcmp(opt.linear_solver,'backslash')
        str = [str sprintf('alpha = C0'' \\ b;\n')];
    elseif strcmp(opt.linear_solver,'qr')
        str = [str sprintf('[Q,R,e] = qr(C0'',0);\n')];
        str = [str sprintf('r = abs(diag(R));\n')];
        str = [str sprintf('tol = max(size(C0)) * eps(r(1));\n')];
        str = [str sprintf('rk = sum(r > tol);\n')];
        str = [str sprintf('alpha = R(1:rk,:) \\ (Q(:,1:rk)''*b);\n')];
        str = [str sprintf('alpha(e,:) = alpha;\n')];
    elseif strcmp(opt.linear_solver,'pinv')
        str = [str sprintf('alpha = pinv(C0'')*b;\n')];
    end
    str = [str sprintf('RR = [alpha''*C1;eye(%d)];\n',length(template.basis))];

end

str = [str cg_format_vector('AM_ind',template.AM_ind,opt.max_str_len)];
str = [str sprintf('AM = RR(AM_ind,:);\n')];

end

function str = generalized_eigenvalue_solver(solv,template,opt)
% Experimental stuff.

str = sprintf('C00 = C0(:,1:end-%d);\nC01 = C0(:,end-%d:end);\n\n',length(template.reducible),length(template.reducible)-1);

str = [str sprintf('[Q,R,e] = qr(C00);\nrr = abs(diag(R));\ntol = max(size(C00)) * eps(rr(1));\nrk = sum(rr > tol);\n')];
str = [str sprintf('A0 = Q(:,rk+1:end)''*C1;\nB0s = Q(:,rk+1:end)''*C01;\n')];

ii = find_mon_indices(template.basis, template.action*template.basis);
red_ind = find_mon_indices(template.action*template.basis,template.reducible);
str = [str cg_format_vector('red_ind',red_ind)];
str = [str sprintf('B0 = zeros(size(A0)); B0(:,red_ind) = -B0s;\n')];

A1 = []; B1 = [];
for k = 1:length(ii)
    if ii(k) > 0
        ak = zeros(1,length(template.basis));
        bk = zeros(1,length(template.basis));
        ak(ii(k)) = 1;        
        bk(k) = 1;
        A1 = [A1; ak];
        B1 = [B1; bk];
    end
end
%norm(A0*V0-B0*V0*D0)

sz = size(A1);
str = [str cg_format_vector('A1_ii',find(A1))];
str = [str cg_format_vector('B1_ii',find(B1))];
str = [str sprintf('A1 = zeros(%d,%d); A1(A1_ii) = 1;\n',size(A1,1),size(A1,2))];
str = [str sprintf('B1 = zeros(%d,%d); B1(B1_ii) = 1;\n',size(B1,1),size(B1,2))];
str = [str sprintf('A = [A0;A1]; B = [B0;B1];\n')];
str = [str sprintf('if size(A,1) > %d\n',length(template.basis))];
str = [str sprintf(' [Q,R,e] = qr([A B],0);\n')];
str = [str sprintf(' A = Q(:,1:%d)''*A;\n B = Q(:,1:%d)''*B;\n',length(template.basis),length(template.basis))];
str = [str sprintf('end\n')];

end

