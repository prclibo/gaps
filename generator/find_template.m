function [ A, opt ] = find_template( eqs, solv, opt, id )

if nargin < 4 || isempty(id)
    id = num2str(feature('getpid'));
elseif ~ischar(id)
    id = num2str(id);
end

if ~exist('./cache','dir'), mkdir('cache'); end
fname = ['./cache/find_template.' solv.name '.' id '.m2'];

disp(fname);

basis = solv.basis;
reducible = solv.reducible;

% if exist(fname,'file'), delete(fname); end

fname_A = ['cache/matrix.' id '.txt'];
fid = fopen(fname,'w');

unk_vars = solv.vars;
header = generate_ring_header(numel(unk_vars), opt );
fprintf(fid,header);

fprintf(fid, 'red = matrix(R, {{%s}});\n', strjoin(string(reducible, false), ', '));
fprintf(fid,'b = matrix(R, {{%s}});\n',strjoin(string(basis, false), ', '));

for i=1:length(eqs)
    fprintf(fid,'eq%u = %s\n',i,string(eqs(i), false));
end
eqname = strjoin(strcat('eq', string(1:numel(eqs))), ', ');
fprintf(fid,'eqs = matrix(R, {{%s}});\n', eqname);
fprintf(fid,'I = ideal eqs;\n');
fprintf(fid,'gbTrace = %d;\n',opt.M2_gbTrace);    

if ~isempty(opt.saturate_mon)
    fprintf(fid,'I0 = I;\n');
    fprintf(fid,'satmon = %s;\n',char(opt.saturate_mon,0));
    fprintf(fid,'I = saturate(I0,satmon);\n');
end

% Find normal set
fprintf(fid,'Q = R/I;\nb0 = lift(basis Q,R);\nuse R\n');

% Express basis in normal set
fprintf(fid,'S = (coefficients(b%%I,Monomials => b0))_1;\n');

fprintf(fid,'if numcols b0 <= numcols b then (\n');
fprintf(fid,'  Sinv = transpose(S)*inverse(S*transpose(S));\n');
fprintf(fid,') else (\n');
fprintf(fid,'  Sinv = inverse(transpose(S)*S)*transpose(S);\n');
fprintf(fid,')\n');

% Construct action matrix
fprintf(fid,'AM = Sinv*((coefficients(red%%I,Monomials => b0))_1);\n');

% Construct target polynomials
fprintf(fid,'pp = red - b*AM;\n');

if ~isempty(opt.saturate_mon)  
    % Find N which lifts target polynomials into the original ideal    
    fprintf(fid,'degs = matrix({apply((entries pp)_0,p->(N=0;while(satmon^N*p %% I0 != 0) do (N=N+1);N))});\n');
    fprintf(fid,'pp = matrix({toList apply(0..numcols pp-1, i->satmon^(degs_(0,i))*pp_(0,i))});\n');
    fname_sat = [ 'cache/saturate_degree.' id '.txt'];        
    fprintf(fid,'%s\n',[' "' fname_sat '" << toString degs << close;']);
end

% Express target polynomials in the generators
fprintf(fid,'A = pp // eqs;\n');

fprintf(fid,'gbRemove(I);\n');


if opt.syzygy_reduction
    fprintf(fid,'M = kernel eqs;\n');
    fprintf(fid,'A = A %% M;\n');        
end

fprintf(fid,'%s\n',[' "' fname_A '" << toString A << close;']);

fprintf(fid,'quit();\n');
fclose(fid);

code = system([opt.M2_path ' ' fname ' --stop']);

if code ~= 0
    warning('M2 return code is %d!', code);
    A = repmat({[]}, numel(eqs), numel(solv.reducible));
    return
end


% read result
if opt.fast_monomial_extraction
    A = extract_monomials(fname_A,length(eqs),length(reducible),numel(unk_vars));
else
    AA = read_m2_matrix(fname_A, m2_vars, unk_vars);
    A = cell(size(AA));
    for k = 1:numel(A)
        [C,mon] = polynomials2matrix(AA(k));
        if C ~= 0
            A{k} = monvec2matrix(mon);
        else
            A{k} = zeros(nvars(eqs(1)),0);
        end
    end
end

% if exist(fname,'file')
%     delete(fname)
% end


if ~isempty(opt.saturate_mon) && exist(fname_sat,'file')    
    opt.saturate_degree = readM2matrix(fname_sat,0);        
    delete(fname_sat);
end






end
