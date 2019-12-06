function [ scheme, norm_scheme ] = find_extraction_scheme( solv, template, opt )

% Figure out which monomials we will have available

sources = {};
% add action
sources{end+1}.mon = template.action;
sources{end}.type = 'action';
sources{end}.idx = 0;

% Add the basis mons
for k = 1:length(template.basis)
    sources{end+1}.mon = template.basis(k);
    sources{end}.type = 'basis';
    sources{end}.idx = k;
end

% Add the basis mons
for k = 1:length(template.reducible)
    sources{end+1}.mon = template.reducible(k);
    sources{end}.type = 'reducible';
    sources{end}.idx = k;
end
    

% Start by finding eigen normalization scheme
norm_scheme = find_eigen_normalization(sources);

% Find variable extraction scheme
scheme = extract_variables(sources);


end


function scheme = extract_variables(sources)
% gather monomials available from start
available_mons = cellfun(@(x) x.mon,sources,'UniformOutput',0);
available_mons = [available_mons{:}]';

target_mons = create_vars(nvars(available_mons(1)));

% remove 1 from sources
ind = find(sum(abs(monvec2matrix(available_mons)))==0);
if ~isempty(ind)
    available_mons(ind) = [];
    sources(ind) = [];
end

% Setup
scheme = cell(0,1);
n_target = length(target_mons);
targets = monvec2matrix(target_mons);
found_target = false(1,length(target_mons));
target_sources = cell(1,n_target);
for k = 1:n_target
    target_sources{k}.mon = target_mons(k);
    target_sources{k}.type = 'target';
    target_sources{k}.idx = k;
end

% Try to express each of the target monomials
done = false;
while ~done
    avail = monvec2matrix([target_mons(found_target);available_mons]);
    avail_src = [target_sources(found_target) sources];
    
    n_avail = size(avail,2);
    
    done = true;
    for tk = find(~found_target)
        target = targets(:,tk);
        
        
        %  1. Check if already available
        [d,i] = min(sum(abs(avail-target*ones(1,n_avail)),1));
        if d == 0
            % found target
            scheme{end+1}.method = 'direct';
            scheme{end}.target_idx = tk;
            scheme{end}.target = target_mons(tk);
            scheme{end}.sources = avail_src(i);
%             fprintf('target=%s exist in available monomials.\n',char(target_mons(tk)));
            found_target(tk) = 1;
            done = false;
            break;
        end
        
        % 2. Try to extract it as a division
        cand_num = find(all(avail-target*ones(1,n_avail) >= 0));
        ind = [];
        for k = cand_num
            % try to find reduction from tmp to zero
            ind = reduce_to_zero(avail(:,k)-target,avail);
            
            if ~isempty(ind)
                % found target
                found_target(tk) = 1;
                
                
                scheme{end+1}.method = 'division';
                scheme{end}.target_idx = tk;
                scheme{end}.target = target_mons(tk);
                scheme{end}.sources = avail_src([k ind]);
                
%                 fprintf('target=%s exist as %s / ',char(target_mons(tk)),char(multipol(1,avail(:,k))));
%                 for ii = ind
%                     fprintf('%s ',char(multipol(1,avail(:,ii))));
%                 end
%                 fprintf('\n');
                done = false;
                break;
            end
            
        end
        if ~isempty(ind)
            break;
        end
        
        % 3. Try to extract it as a sqrt
        target2 = target*2;
        [d,i] = min(sum(abs(avail-target2*ones(1,n_avail))));
        if d == 0
            % found target
            scheme{end+1}.method = 'sqrt';
            scheme{end}.target_idx = tk;
            scheme{end}.target = target_mons(tk);
            scheme{end}.sources = avail_src(i);
%             fprintf('target=%s exist as sqrt in avail.\n',char(target_mons(tk)));
            found_target(tk) = 1;
            done = false;
            break;
        end
    end
end

if any(~found_target)
%     str = 'unable to extract: ';
%     for k = find(~found_target)
%         str = [str sprintf('%s ',char(target_mons(k)))];
%     end
%     warning(str);
end



end

function [ norm_scheme ] = find_eigen_normalization( sources )

norm_scheme = [];

% Check if 1 is available
available_mons = cellfun(@(x) x.mon,sources,'UniformOutput',0);
available_mons = [available_mons{:}]';
avail = monvec2matrix(available_mons);
[d,i] = min(sum(abs(avail)));
if d == 0
    norm_scheme.method = 'direct';
    norm_scheme.sources = sources{i};
    return;
end

% locate action mon
ind = [];
for k = 1:length(sources);
    if strcmp(sources{k}.type,'action')
        ind = k;
        break;
    end
end
if isempty(ind)
    error('No action monomial?');
end

action = sources{ind};
available_mons = cellfun(@(x) x.mon,sources,'UniformOutput',0);
available_mons = [available_mons{:}]';

% remove action from sources
available_mons(ind) = [];
sources(ind) = [];
avail = monvec2matrix(available_mons);

% Try to scale with eigenvalue
tt = monvec2matrix(action.mon);
ind = reduce_to_zero(tt,avail);

if ~isempty(ind)
    % Success!
    norm_scheme.method = 'eigen';
    norm_scheme.sources = sources(ind);
    return;
end

ind = reduce_to_zero(tt/2,avail);
if ~isempty(ind)
    % Success!
    norm_scheme.method = 'sqrt_eigen';
    norm_scheme.sources = sources(ind);
    return;
end

end




function ind = reduce_to_zero(tt,avail)
% find candidates
divisible = find(all(tt*ones(1,size(avail,2))-avail>=0));
if isempty(divisible)
    ind = [];
    return;
end

for k = divisible
    % check if we are done
    if sum(abs(tt-avail(:,k))) == 0
        ind = k;
        return;
    end
    
    ind = reduce_to_zero(tt-avail(:,k),avail);
    if ~isempty(ind)
        ind = [k ind];
        return;
    end
end

end