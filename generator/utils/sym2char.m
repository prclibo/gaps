function c = sym2char(s, paren)

c = char(s);
if numel(s) > 1
    c = c(8:end - 1);
end

% Replace '], [' by ']; ['
if nargin == 1 || strcmp(paren, '[]')
    c = strrep(c, '], [', ']; [');
end

if nargin > 1
    if isempty(paren)
        c = strrep(c, '[', '');
        c = strrep(c, ']', '');
    else
        assert(any(ismember(paren, {'{}', '()'})));
        c = strrep(c, '[', paren(1));
        c = strrep(c, ']', paren(2));
    end
end