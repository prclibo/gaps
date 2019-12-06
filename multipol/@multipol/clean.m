function p1 = clean(p1);
% MULTIPOL/CLEAN operator
% Remove monomials with coefficent zero
% p1 = clean(p1);


    
I = find(p1.coeffs==0); 
p1.coeffs(I)=[];
p1.monomials(:,I)=[];
