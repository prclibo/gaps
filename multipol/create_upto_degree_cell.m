function p1 = create_upto_degree_cell(vars_degree,max_degree)
% p1 = create_upto_degree(vars_degree,max_degree)
% CREATE_UPTO_DEGREE - Creates a multipol cell object with degree up to
% the numbers in var degree. The maximum degree accepted is in
% max_degree. 
%   

if nargin == 1
   max_degree = inf;
end

mons = [];
for ii = 1:length(vars_degree)
   mons = [mons 0:vars_degree(ii)];
end
mons = nchoosek(mons',length(vars_degree));

for ii = 1:length(vars_degree)
   to_big = find(mons(:,ii)>vars_degree(ii));
   mons(to_big,:) = [];
end

to_big = find(sum(mons')>max_degree);
mons(to_big,:) = [];

coeffs = ones(1,size(mons,1));
p1=multipol(coeffs,mons');
p1=sort(p1);

m = monomials(p1);
clear p1;
for k = 1 : size(m, 2)
    p1{k} = multipol(1, m(:, k));
end
