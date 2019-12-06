function p3 = subs(p1,var,p2);
% MULTIPOL/SUBS operator
% Substitution
% p3 = subs(p1,var,p2);
% p3 = p1(..,p2,..); where p2 is one place var.
% p2 is a double or a multipol object.

iscell = false;
if isa(p1, 'cell')
    iscell = true;
    for k = 1 : length(p1);
        polys(k) = p1{k};
    end
else
    polys = p1;
end

for k = 1 : length(polys)
    p1 = polys(k);
    if isa(p2,'double');
        coeffs=p1.coeffs;
        monomials=p1.monomials;
        places = find(p1.monomials(var,:));
        coeffs(places)=p1.coeffs(places).*(p2).^p1.monomials(var,places);
        monomials(var,:) = [];

        p3(k) = multipol(coeffs,monomials);
        p3(k) = sort(p3(k));
        p3(k) = squeeze(p3(k));

    else
        powers = p1.monomials(var,:);
        p1.monomials(var,:) = 0;
        p3(k) = multipol;
        for jj=1:size(p1.monomials,2)
            p3(k) = p3(k) + multipol(p1.coeffs(jj),p1.monomials(:,jj))...
                * p2^powers(jj);
        end
        p3(k) = squeeze(p3(k));
    end
end

if(iscell)
    pp = p3;
    p3 = cell(length(pp), 1);
    for k = 1 : length(p3)
        p3{k} = pp(k);
    end
end