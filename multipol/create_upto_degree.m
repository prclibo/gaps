function p1 = create_upto_degree(vars_degree,max_degree)
% p1 = create_upto_degree(vars_degree,max_degree)
% CREATE_UPTO_DEGREE - Creates a multipol object with degree up to
% the numbers in var degree. The maximum degree accepted is in
% max_degree. 
%   

if nargin == 1
   max_degree = inf;
end

monomials = [];
for ii = 1:length(vars_degree)
   monomials = [monomials 0:vars_degree(ii)];
end
monomials = nchoosek(monomials',length(vars_degree));

for ii = 1:length(vars_degree)
   to_big = find(monomials(:,ii)>vars_degree(ii));
   monomials(to_big,:) = [];
end

to_big = find(sum(monomials')>max_degree);
monomials(to_big,:) = [];

coeffs = ones(1,size(monomials,1));
p1=multipol(coeffs,monomials');
p1=sort(p1);
