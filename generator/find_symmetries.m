function [ cc, pp ] = find_symmetries( eqs, max_p )
% Automatically finds any axis aligned symmetries
% TODO: Make this less of a hack
[C,mon] = polynomials2matrix(eqs);
monv = monvec2matrix(mon);
ind = abs(C)>0;
cc = [];
pp = [];

for p = 2:max_p
    c = next(zeros(size(monv,1),1),p);
    while ~isempty(c)
        
        good = true;

        % check for symmetry
        for k = 1:length(eqs)
            if all(~ind(k,:))
                % zero equation, continue
                continue;
            end
            
            rem = mod(c'*monv(:,ind(k,:)),p);
            if any(rem~=rem(1))
                good = false;
                break;
            end
        end
        
        if good
            % found symmetry
            cc = [cc c];
            pp = [pp p];
        end
        c = next(c,p);
    end
end

% try to remove factor symmetries
%  i.e. (2,0) p=4 = (1,0) p=2
good = true(1,length(pp));
for k = 1:length(pp) 
    ck = cc(:,k);
    pk = pp(k);
    for fk = factor(pk);
        if norm(floor(ck/fk)*fk-ck) == 0
            good(k) = 0;
        end
    end
end
cc = cc(:,good);
pp = pp(:,good);


end

function c = next(c,p)
    c(1) = c(1)+1;
    k = 1;
    while c(k) == p
        c(k) = 0;
        k = k + 1;
        if k > length(c)
            c = [];
            return;
        end
        c(k) = c(k)+1;
    end
end

