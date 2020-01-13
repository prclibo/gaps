function [ result ] = setup_reduced_eigenvector_solver( template )

action = template.action;
basis = template.basis;

% Setup dummy AM
RR = [ones(length(template.reducible),length(basis)); eye(length(basis))];
AM = RR(template.AM_ind,:);

eqs = action*basis(:)-AM*basis(:);
ind_non_trivial = find(template.AM_ind <= length(template.reducible));

nv = max(nvars(basis));

av = monomials(action);
assert(sum(av)==1, 'NYI: reduced_eigenvector_solver currently only works for single variable actions.');
av = find(av);

my_vars = vars(basis(1));
xx = multipol(my_vars(:), 'vars', my_vars);
xx(av) = multipol(1, 'vars', my_vars);

% Reduced basis (i.e. after subst. action)
b0 = monvec(sum(evaluate(basis,xx)));
b0v = monvec2matrix(b0);

% Make sure we have 1 in b0
assert(any(sum(b0v)==0),'NYI: reduced_eigenvector_solver currently only works for cases where 1 is in the reduced basis.')

% check that we have all monomials that we need
ind_var = zeros(nv,1);
ind_var(av) = -1;
for k = 1:length(b0)
    if sum(b0v(:,k))==1
        ind_var(find(b0v(:,k))) = k;
    end
end
assert(all(ind_var~=0),'NYI reduced_eigenvector_solver currently only handles cases where are variables are trivially available.')

% Compute monomials in reduced eigenvector problem
A = collect_terms(basis,b0);

result = [];
result.AA_sz = [length(ind_non_trivial),length(b0)];
result.b0 = b0;
result.ind_var = ind_var;
result.ind_unit = find(sum(b0v)==0);
assert(result.ind_unit == result.AA_sz(2),'NYI: reduced_eigenvector_solver currently only handles cases where 1 is last in basis.')
result.ind_non_trivial = ind_non_trivial;
result.ind_AA1 = {};
result.ind_AA2 = [];
result.max_power = nan;

assert(result.AA_sz(1) >= result.AA_sz(2)-1,'reduced_eigenvector_solver something is strange?')

for k = 1:size(A,2)
    ii = find(A(:,k)~=0);
    if isempty(ii)
        continue;
    end
%     str = '';
%     for i = ii'
%         if A(i,k) == 1
%             str = [str sprintf('AM(:,%d)+',i)];
%         else
%             str = [str sprintf('%s*AM(:,%d)+',char(A(i,k)),i)];
%         end
%     end
%     fprintf('AA(:,%d) = %s;\n',k,str(1:end-1));
    result.ind_AA1{k} = [];
    for i = ii'
        mv = monomials(A(i,k));
        result.ind_AA1{k}(:,end+1) = [i;sum(mv)];            
    end
end
  
for k = 1:size(A,2)
 ii = find(A(:,k)~=0);
 
    if isempty(ii)
        continue;
    end

%     str = '';
    for i = ii'
        if ~ismember(i,ind_non_trivial)
            continue;
        end        
        j = find(ind_non_trivial == i);
        
        mv = monomials(action*A(i,k));
        result.ind_AA2(:,end+1) = [j;sum(mv)];
%         if A(i,k) == 1
%             str = [str sprintf('AA(%d,%d) = AA(%d,%d) - y;\n',j,k,j,k)];
%         else
%             str = [str sprintf('AA(%d,%d) = AA(%d,%d) - %s;\n',j,k,j,k,char(action*A(i,k)))];
%         end
    end
%     fprintf(str);
end

result.max_power = max(result.ind_AA2(2,:));

end

