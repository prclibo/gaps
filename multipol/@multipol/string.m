function s = string(mp,lex_vars)
% MULTIPOL/STRING
% Convert polynomial expression to readable string representation
if nargin<2
    lex_vars = true;
end

s = repmat("", size(mp));

for k = 1:numel(mp)
    if lex_vars
        unk_vars = mp(k).vars(:);
    else
        n = numel(mp(k).vars);
        unk_vars = sym('x%d', [n, 1]);
    end
    
    monos = string(prod(unk_vars .^ int8(mp(k).monomials)));
    one_mask = logical(mp(k).coeffs == 1);
    coefs = string(mp(k).coeffs);
    coefs(one_mask) = "";
    mul_symb = repmat("*", size(coefs));
    mul_symb(one_mask) = "";
    terms = strcat(coefs, mul_symb, monos);
    line = strjoin(terms, ' + ');
    s(k) = line;
end