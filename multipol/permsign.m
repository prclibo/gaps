function s = permsign(p)
% s = permsign(p) returns the sign of the permutation p
s = 1;
for i = 1 : length(p)-1
    for j = i+1 : length(p)
      if(p(i)>p(j))
          s = -s;
      end
    end
end