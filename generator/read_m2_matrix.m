function [ A ] = read_m2_matrix( fname, unk_vars )


fid = fopen(fname,'r');
line = fgetl(fid);
line = strrep(line,'matrix ','');
line = strrep(line,'{','[');
line = strrep(line,'}','];');
fclose(fid);

if nargin > 1
    for i = 1:numel(unk_vars)
        cmd = sprintf('x%d = multipol(unk_vars(i), ''vars'', unk_vars);', i);
        eval(cmd);
    end
end
A = eval(line);
if nargin > 1
    A = A * multipol(1, 'vars', unk_vars);
end
end

