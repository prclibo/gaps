function vals = subs_var(eqs, var_subs, varargin)

rng(21);
verbose_mask = strcmp(varargin, 'verbose');
verbose = any(verbose_mask);
varargin(verbose_mask) = [];

parser = inputParser;
parser.KeepUnmatched = true;
parser.addOptional('zp', 0);
parser.parse(varargin{:});

zp_p = parser.Results.zp;
is_zp = zp_p > 0;

if isa(eqs, 'sym')
    eqs_ = {eqs};
elseif isstruct(eqs)
    eqs_ = struct2cell(eqs);
else
    eqs_ = eqs;
    assert(iscell(eqs), 'Wrong type of eqs');
end

var_names = fieldnames(var_subs);
subs_ = struct2cell(var_subs);
for i = 1:numel(var_names)
    if is_zp
        assign(var_names{i}, subs(subs_{i}));
    else
        assign(var_names{i}, subs(subs_{i}));
    end
end

subs_numel = cellfun(@(x) numel(x), subs_);
if any(subs_numel > 1) && ~iscell(eqs)
    error('Found non-scalar subs. You should evaluate cell as the output shape varies.');
end
% kwn_vars = sym(var_names);
% all_vars = symvar([eqs_{:}]);
% unk_vars = setdiff(all_vars, kwn_vars);

vals = cell(size(eqs_));
for ci = 1:numel(eqs_)
    if verbose
        fprintf('...In total %d terms\n', numel(eqs));
        fprintf('...Evaluating instances for           ');
    end
    
    if is_zp, vals{ci} = sym(zeros(size(eqs_{ci})));
    else, vals{ci} = zeros(size(eqs_{ci})); end

    vals{ci} = subs(eqs_{ci});
    
    if isempty(symvar(vals{ci}))
        if is_zp
            vals{ci} = mod(vals{ci}, zp_p);
        else
            vals{ci} = double(vals{ci});
        end
    else
        indices = find(logical(eqs_{ci} ~= sym(0)));
        for i = indices(:)'
            if verbose; fprintf('\b\b\b\b\b\b\b\b%8d', i); end
            [C, T] = coeffs(vals{ci}(i));
            if is_zp
                C = mod(C, zp_p);
            else
                C = double(C);
            end
            vals{ci}(i) = dot(C, T);
        end
        
    end
    
    if verbose; fprintf('\n'); end
end

if isstruct(eqs)
    vals = cell2struct(vals, fieldnames(eqs));
elseif isa(eqs, 'sym')
    vals = vals{1};
end

end

function assign(name, val)
assignin('caller', name, val);
end
