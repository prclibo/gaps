function [in_zp, out_zp] = vzp_rand_relpose_pinhole(N, p, in_fet, out_fet, varargin)

parser = inputParser;
parser.addParameter('ZeroScrewTransl', 'none',...
    @(x) any(validatestring(x, {'inner', 'outer', 'none'})));
parser.parse(varargin{:});

quat = zp_rand_unit(4, p);
s = quat(1);
u1 = quat(2);
u2 = quat(3);
u3 = quat(4);
u = [u1; u2; u3];

R = zp_quat2dcm(quat, p);

tau = 4 * s^2 - 1; tau = mod(tau, p);

v1 = u1 * zp_inv(s, p); v1 = mod(v1, p);
v2 = u2 * zp_inv(s, p); v2 = mod(v2, p);
v3 = u3 * zp_inv(s, p); v3 = mod(v3, p);
v = [v1; v2; v3];

t = zp_rand_unit(3, p);

if any(strcmpi(parser.Results.ZeroScrewTransl, {'inner', 'outer'}))
    error('Not implemented yet');
    t = skew(u) * t; t = mod(t, p);
end

Q = sym(randi(p, [3, N]));
E = skew(t) * R;
QQ = R * Q + repmat(t, [1, N]);
q = Q; qq = QQ;

if any(ismember('NE', in_fet))
    AE = sym(zeros([N, 9]));
    for i = 1:N
        AE(i, :) = reshape(qq(:, i) * q(:, i)', 1, []);
    end
    if strcmpi(parser.Results.ZeroScrewTransl, 'inner')
        AE(N + 1, :) = [1, 0, 0, 0, 1, 0, 0, 0, 1];
    end
    
    NE = null(AE);
    w = NE \ E(:);
    
    NE = zp_eval_rational(NE, p);
    
    w = w .* zp_inv(w(end), p); w = mod(w, p);
    w = w(1:end - 1);
    
    NE = reshape(NE, 3, 3, size(NE, 2));
    NE = permute(NE, [3, 1, 2]);
    E = mod(E, p);
end
q = mod(q, p);
qq = mod(qq, p);

in_zp = pack_pars(in_fet);
out_zp = pack_pars(out_fet);

end

function ret = pack_pars(fetches)
ret = struct();
for var_name = fetches
    ret.(var_name{1}) = evalin('caller', var_name{1});
end

end
% clear varargin i parser AE QQ Q;
% var_names = who();
% var_zp = repmat(struct('var', [], 'subs', []), size(var_names));
% for i = 1:numel(var_names)
%     name = var_names{i};
%     var_zp(i).var = sym(name);
%     var_zp(i).subs = eval(name);
% end
