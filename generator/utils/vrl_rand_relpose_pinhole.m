function [in_rl, out_rl] = vrl_rand_relpose_pinhole(N, in_fet, out_fet, varargin)

parser = inputParser;
parser.addParameter('ZeroScrewTransl', 'none',...
    @(x) any(validatestring(x, {'inner', 'outer', 'none'})));
parser.parse(varargin{:});

quat = normc(rand([4, 1]));
s = quat(1);
u1 = quat(2);
u2 = quat(3);
u3 = quat(4);
u = [u1; u2; u3];

R = sym_quat2dcm(quat);

tau = 4 * s^2 - 1;

v1 = u1 / s;
v2 = u2 / s;
v3 = u3 / s;
v = [v1; v2; v3];

t = normc(rand([3, 1]));

if any(strcmpi(parser.Results.ZeroScrewTransl, {'inner', 'outer'}))
    error('Not implemented yet');
    t = skew(u) * t;
end

Q = rand([3, N]);
E = skew(t) * R;
QQ = R * Q + repmat(t, [1, N]);
q = Q; qq = QQ;

if any(ismember('NE', in_fet))
    AE = zeros([N, 9]);
    for i = 1:N
        AE(i, :) = reshape(qq(:, i) * q(:, i)', 1, []);
    end
    if strcmpi(parser.Results.ZeroScrewTransl, 'inner')
        AE(N + 1, :) = [1, 0, 0, 0, 1, 0, 0, 0, 1];
    end
    
    NE = null(AE);
    w = NE \ E(:);
        
    w = w ./ w(end);
    w = w(1:end - 1);
    
    NE = reshape(NE, 3, 3, size(NE, 2));
    NE = permute(NE, [3, 1, 2]);
end

in_rl = pack_pars(in_fet);
out_rl = pack_pars(out_fet);

end

function ret = pack_pars(fetches)
ret = struct();
for var_name = fetches
    ret.(var_name{1}) = evalin('caller', var_name{1});
end

end
% clear varargin i parser AE QQ Q;
% var_names = who();
% var_rl = repmat(struct('var', [], 'subs', []), size(var_names));
% for i = 1:numel(var_names)
%     name = var_names{i};
%     var_rl(i).var = sym(name);
%     var_rl(i).subs = eval(name);
% end
