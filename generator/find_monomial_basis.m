function [ basis ] = find_monomial_basis(eqs_zp, solv, opt, id )

if nargin < 3 || isempty(id)
    id = num2str(feature('getpid'));
elseif ~ischar(id)
    id = num2str(id);
end

if ~exist('./cache','dir'), mkdir('cache'); end

fname = ['./cache/find_monomial_basis.' solv.name '.' id '.m2'];
disp(fname);
if exist(fname,'file'), delete(fname); end

fname_b = ['cache/basis.' id '.txt'];
fid = fopen(fname,'w');

header = generate_ring_header(nvars(eqs_zp(1)), opt );
fprintf(fid,header);

for i=1:length(eqs_zp)
    fprintf(fid,'eq%u = %s\n',i,string(eqs_zp(i), false));
end
eqname = strjoin(strcat('eq', string(1:numel(eqs_zp))), ', ');
fprintf(fid,'eqs = {%s}\n',eqname);
fprintf(fid,'I = ideal eqs;\n');
fprintf(fid,'gbTrace = %d;\n',opt.M2_gbTrace);

if isempty(opt.saturate_mon)
    fprintf(fid,'Q = R/I;\n');
    fprintf(fid,'b = basis Q;\n');
else
    fprintf(fid,'satmon = %s\n',char(opt.saturate_mon,0));
    fprintf(fid,'J = saturate(I,satmon);\n');
    fprintf(fid,'Q = R/J;\n');
    fprintf(fid,'b = basis Q;\n');
end

fprintf(fid,'%s\n',[' "' fname_b '" << toString b << close']);
fprintf(fid,'quit();\n');
fclose(fid);

code = system([opt.M2_path ' ' fname ' --stop']);

if code ~= 0
    warning('M2 return code is %d!', code);
    basis = [];
    return
end

% read result
basis = read_m2_matrix(fname_b, vars(eqs_zp(1)));
assert(~isempty(basis), 'Didnot find basis! Check if the ideal is zero-dimensional');

% if exist(fname,'file'), delete(fname); end
% if exist(fname_b,'file'), delete(fname_b); end