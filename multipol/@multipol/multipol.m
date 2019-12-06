function s = multipol(x, varargin)

if nargin==0
	s.coeffs = 0;
	s.monomials = 0;
    s.vars = sym([]);
    s.str = "";
    s = class(s,'multipol');
else
    mono_degs = [];
    if ~isempty(varargin) && ~ischar(varargin{1})
        mono_degs = varargin{1};
        varargin(1) = [];
    end
    parser = inputParser;
    parser.addParameter('vars', sym([]), @(x) isa(x, 'sym') || isa(x, 'multipol'));
    parser.parse(varargin{:});
    vars = parser.Results.vars;
    if isa(vars, 'multipol'), vars = vars.vars; end
    if isempty(x), x = 0; end
    assert(isempty(mono_degs) || all(numel(x) == size(mono_degs, 2)),...
        'Size not match!');
    if isa(x, 'multipol')
        if numel(x)==1
            s.coeffs = x.coeffs;
            s.monomials = x.monomials;
            s.vars = x.vars;
            s.str = x.str;
            if ~isempty(vars)
                s.vars = vars;
            end
            s = class(s,'multipol');
            if ~isempty(vars)
                s.str = x.str;
            end
        else
            % Should do test. What is the scenario for this?
            for i=numel(x):-1:1
                if isempty(vars)
                    s(i) = multipol(x(i));
                else
                    s(i) = multipol(x(i), 'vars', vars);
                end
            end
            s = reshape(s,size(x));
        end
    elseif isa(x, 'sym') && isempty(mono_degs)
        if isempty(vars)
            vars = symvar(x);
        end
        s = from_sym(x, vars);
%         s = class(s,'multipol');
    else
        assert(isnumeric(x) || isa(x, 'sym'));
        if ~isempty(mono_degs)
            if isempty(vars)
                vars = sym('x%d', [size(mono_degs, 1), 1]);
            else
                assert(numel(vars) == size(mono_degs, 1), 'Sizes are %d, %d', numel(vars), size(mono_degs, 1));
            end
            if isnumeric(x)
                coeff_vars = [];
            else
                coeff_vars = symvar(x);
            end
            assert(isempty(coeff_vars) || isempty(intersect(coeff_vars, vars)),...
                'Coeff should not contain poly vars.');
        else
            mono_degs = zeros(numel(vars), 1);
        end
        s = repmat(multipol, size(x, 1), 1);
        for i = 1:numel(s)
            coefs = x(1, :);
            nz_mask = coefs ~= 0;
            s(i).coeffs = coefs(nz_mask);
            s(i).monomials = mono_degs(:, nz_mask);
            s(i).vars = vars;
            s(i).str = string(s(i));
        end
    end
end
end
%%

function s = from_sym(eqs, vars)

s = repmat(multipol, size(eqs));
var_str = evalc('disp(vars)');
for i = 1:numel(eqs)
    eq = eqs(i);
    % Same speed the below.
    % pl = mupadmex('poly2list', char(eq), var_str);
    pl = feval(symengine, 'poly2list', eq, var_str);
    nterms = numel(pl);
    m = zeros(numel(vars),nterms);
    cidxs = zeros(1,nterms);
    coefs = sym(zeros(1,nterms));
    pl = feval(symengine, 'matrix', pl);
    coefs = pl(:, 1).';
    for j=1:nterms
        m(:,j) = pl(j, 2);
    end
    s(i) = multipol(coefs, m, 'vars', vars);
end

% s = repmat(multipol, size(eqs));
% var_str = evalc('disp(vars)');
% for i = 1:numel(eqs)
%     eq = eqs(i);
%     
%     [coefs, T] = coeffs(eq, vars);
%     nterms = numel(T);
%     m = zeros(numel(vars),nterms);
%     for j = 1:nterms
%         pl = mupadmex('poly2list', char(T(j)), var_str);
%         pl = pl(1);
%         m(:, j) = pl(2);
%     end
%     
%     s(i) = multipol(coefs,m, 'vars', vars);
% end

% s = repmat(multipol, size(eqs));
% var_str = evalc('disp(vars)');
% for i = 1:numel(eqs)
%     eq = eqs(i);
%     
%     [coefs, T] = coeffs(eq, vars);
%     nterms = numel(T);
%     m = zeros(numel(vars),nterms);
%     for j = 1:nterms
% %         pl = mupadmex('poly2list', char(T(j)), var_str);
%         pl = feval(symengine, 'degreevec', T(j), var_str);
%         m(:, j) = pl;
%     end
%     
%     s(i) = multipol(coefs,m, 'vars', vars);
% end

end
