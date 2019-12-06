function f = evaluate(p1,x);
% f = evaluate(p1,x);
% MULTIPOL/EVALUATE operator
% Calculates the value of the polynomial p1
% at the point x

if(numel(p1)>1)
    m = size(p1,1);
    n = size(p1,2);
    for i = 1 : m
        for j = 1 : n
            f(i,j) = evaluate(p1(i,j),x);
        end
    end
else
    if(nterms(p1)==0)
        f = 0;
        return;
    end
    c=p1.coeffs;
    m=p1.monomials;
    %keyboard;
    f=c;
    for i=1:size(m,1);
        f=f.*(x(i,1).^m(i,:));
    end;
    f=sum(f);
end