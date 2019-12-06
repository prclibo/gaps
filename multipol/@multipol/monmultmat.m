function monmat = monmultmat(p1,p2);
% MULTIPOL/MONMULTMAT 
% Calculates the matrices that are used to multiplicate with a
% monomial. This function is only dependent of the monomial structure and
% not the coefficents.
% 
%    monmat = monmultmat(p1,p2);
% The degree of p1 should be greater than p2.

monomials1=p1.monomials;
monomials2=p2.monomials;
degree_diff=max(monomials1')-max(monomials2');
monmat = [];
nbr_terms1 = size(monomials1,2);
nbr_terms2 = size(monomials2,2);
degree_diff = flipud(unique(nchoosek([0:degree_diff(1) 0:degree_diff(2)],2),'rows'));

for ii = 1:size(degree_diff,1)
   monmat_tmp = zeros(nbr_terms2,nbr_terms1);
   monomials2_tmp = monomials2+degree_diff(ii,:)'*ones(1,size(monomials2,2));
   monomials2_tmp(:,end+1) = -1e10*ones(size(monomials2,1),1);
   % Do this for all monomials
   nn = 1;
   for kk = 1:nbr_terms1
      if sum(abs(monomials1(:,kk)-monomials2_tmp(:,nn)))==0
	 monmat_tmp(nn,kk) = 1;
	 nn = nn+1;
      end
   end
   monmat{end+1}=monmat_tmp;
end
