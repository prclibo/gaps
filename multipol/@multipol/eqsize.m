function varargout = eqsize(varargin)
% MULTIPOL/EQSIZE Resize internal monomial matrices to consistent size between all input objects

if nargin == 1
    varargout = varargin;
elseif nargin == 2
    if numel(varargin{1}(1).vars) == numel(varargin{2}(1).vars) &&...
            all(varargin{1}(1).vars(:) == varargin{2}(1).vars(:))
        varargout = varargin;
        return;
    end
    unk_vars = union(varargin{1}(1).vars, varargin{2}(1).vars, 'stable');
    varargout = cell(1, 2);
    for i = 1:2
        p = varargin{i};
        varargout{i} = repmat(multipol, size(varargin{i}));
        [~, ia, ib] = intersect(unk_vars, p(1).vars, 'stable');
        for j = 1:numel(p)
            assert(all(p(j).vars(:) == p(1).vars(:)));
            m = zeros(numel(unk_vars), size(p(j).monomials, 2));
            m(ia, :) = p(j).monomials(ib, :);
            varargout{i}(j) = multipol(p(j).coeffs, m, 'vars', unk_vars);
        end
    end
elseif nargin > 2
    tmp = eqsize(varargin{2:end});
    varargout = eqsize(varargin{1}, tmp{1});
end
