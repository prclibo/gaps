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
        unk_vars = sym('x', [n, 1]);
    end
    
    if ~isempty(mp(k).monomials)
        monos = string(prod(unk_vars .^ int8(mp(k).monomials), 1));
        one_mask = logical(mp(k).coeffs == 1);
        coefs = string(mp(k).coeffs);
        coefs(one_mask) = "";
        mul_symb = repmat("*", size(coefs));
        mul_symb(one_mask) = "";
        terms = strcat(coefs, mul_symb, monos);
        line = strjoin(terms, ' + ');
    else
        line = '0';
    end
    s(k) = line;
    

%     monos = sym2charcell(prod(unk_vars .^ int8(mp(k).monomials)));
%     
%     if size(mp(k).coeffs, 2) == 1 && logical(mp(k).coeffs == 1) &&...
%         numel(monos) == 1
%         s(k) = char(monos(1));
%         continue
%     end
%     
%     one_mask = logical(mp(k).coeffs == 1);
%     coefs = sym2charcell(mp(k).coeffs);
%     coefs(one_mask) = {''};
%     mul_symb = repmat({'*'}, size(coefs));
%     mul_symb(one_mask) = {''};
%     terms = strcat(coefs, mul_symb, monos);
%     line = strjoin(terms, ' + ');
%     s(k) = line;
end