function eqs2phcpack(eq,filename)

fid = fopen(filename,'w');
nvar = nvars(eq(1));

fprintf(fid,'%u %u\n',nvar,nvar);

for i=1:numel(eq)
	fprintf(fid,'%s;\n',char(eq(i),0,32));
end
