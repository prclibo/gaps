function eqs2maple(eq,fname)

if nargin<2
	fid = 1;
else
	fid = fopen(fname,'w');
end

fprintf(fid,'with(Groebner):\n');
fprintf(fid,'F:={');
for i=1:numel(eq)-1
	fprintf(fid,'%s,',char(sym(eq(i),0)));
end
fprintf(fid,'%s}:\n',char(sym(eq(end),0)));
fprintf(fid,'m_ord:=tdeg(');
for i=1:nvars(eq(1))-1
	fprintf(fid,'x%u,',i);
end
fprintf(fid,'x%u):\n',nvars(eq(1)));
fprintf(fid,'G:=Basis(F,m_ord):\n');
fprintf(fid,'NormalSet(G,m_ord)[1]\n');

if nargin>1
	fclose(fid);
end