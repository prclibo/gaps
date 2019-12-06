function eqs2m2(eqs,fname)


[C M] = polynomials2matrix(eqs);
C = round(C);
eqs = C*M;

varstr = '';
for i=1:nvars(eqs(1))
	varstr = [varstr sprintf('x%u,',i)]; %#ok
end

fid = fopen(fname,'w');
fprintf(fid,'%s\n','KK = ZZ / 30097;');
fprintf(fid,'%s\n',['R = KK[' varstr(1:end-1) '];']);

eqname = '';
for i=1:length(eqs)
	fprintf(fid,'eq%u = %s\n',i,char(eqs(i),0));
	eqname = [eqname sprintf('eq%u,',i)]; %#ok
end

fprintf(fid,'%s\n',['eqs = {' eqname(1:end-1) '};']);
fprintf(fid,'%s\n','I = ideal eqs;');
fprintf(fid,'%s\n','gbTrace = 3');
fprintf(fid,'%s\n','dim I, degree I');
fprintf(fid,'%s\n','dimi = dim I;');
fprintf(fid,'%s\n','degreei = degree I;');
fprintf(fid,'%s\n','nrsols = degreei');
fprintf(fid,'%s\n','a = gens gb I');
fprintf(fid,'%s\n','lta = leadTerm a');
fprintf(fid,'%s\n',[' "' fname '.out" << dimi << endl << nrsols << endl << lta << endl << close']);
fprintf(fid,'%s\n','quit');
fclose(fid);