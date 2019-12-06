function str = cg_format_vector(name,v,max_length)
if nargin < 3
    max_length = 50;
end
str = sprintf('%s = [',name);
while length(v) > max_length
    str = [str sprintf('%d,',v(1:max_length))];
    str = [str sprintf('...\n')];
    v = v(max_length+1:end);
end
if ~isempty(v)
    str = [str sprintf('%d,',v)];
    str = str(1:end-1);
end
str = [str sprintf('];\n')];
