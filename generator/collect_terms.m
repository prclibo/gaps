function [ A, rem ] = collect_terms( eqs, terms )
[cc,mm] = polynomials2matrix(eqs);

mmv = monvec2matrix(mm);
tv = monvec2matrix(terms);

% ensure that they have the same number of variables
tv = [tv; zeros(max(0,size(mmv,1)-size(tv,1)),size(tv,2))];
mmv = [mmv; zeros(max(0,size(tv,1)-size(mmv,1)),size(mmv,2))];

A = zeros(length(eqs),length(terms));
mzero = multipol(1,zeros(nvars(eqs(1)),1));
A = A * mzero;
rem = zeros(length(eqs),1) * mzero;
nz_rem = false;

for k = 1:length(mm)
    done = false;
    for i = 1:length(terms)
        % check it mm(i) is divisible by tv(i)
        mdiv = mmv(:,k) - tv(:,i);
        if any(mdiv < 0)
            continue;
        end
        
        A(:,i) = A(:,i) + cc(:,k)*multipol(1,mdiv);
        done = true;
        break;
    end
    
    if ~done
        nz_rem = true;
        rem = rem + cc(:,k)*mm(k);
    end
end

if nargout == 1 && nz_rem
    warning('Remainder was nonzero but ignored.')
end


end


