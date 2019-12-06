function f = power(p1,B);
% f = power(p1,b);
% MULTIPOL/POWER operator
% This function raises each element of p1 to the corresponding power in b.
%
%    C = p1.^B
%    C = power(p1,B)

if(numel(p1)>1)
    m = size(p1,1);
    n = size(p1,2);
    if(numel(B)==1)
        for i = 1 : m
            for j = 1 : n
                f(i,j) = power(p1(i,j),B);
            end
        end
    elseif(all(size(B)==size(p1))),
        for i = 1 : m
            for j = 1 : n
                f(i,j) = power(p1(i,j),B(i,j));
            end
        end
    end
    %keyboard;
else
    if(nterms(p1)==0)
        f = 0;
        return;
    end
    if numel(p1)==1,
        if numel(B)>1,
            m = size(B,1);
            n = size(B,2);
            for i = 1 : m
                for j = 1 : n
                    f(i,j) = power(p1,B(i,j));
                end
            end
        elseif numel(B)==1,
            f = p1^B;
        end
    end
end