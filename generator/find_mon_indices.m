function [ ind ] = find_mon_indices( mons, target )

if isa(mons, 'multipol')
    mons = monvec2matrix(mons);
end
if isa(target, 'multipol')
    target = monvec2matrix(target);
end

ind = arrayfun(@(x) find(sum(abs(mons-target(:,x)*ones(1,size(mons,2))),1)==0,1), 1:size(target,2),'UniformOutput',0);
for k = 1:length(ind)
    if isempty(ind{k})
        ind{k} = -1;
    end
end
ind = cell2mat(ind);
end

