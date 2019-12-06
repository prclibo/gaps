function s = eqs2mupad(eq,zp)

if nargin<2

	s = '[';
	for i=1:numel(eq)-1
		s = [s sprintf('%s, ',char(sym(eq(i),0)))];
	end
	s = [s sprintf('%s]',char(sym(eq(end),0)))];
	
else
	
	vars = '[';
	for i=1:nvars(eq(1))-1
		vars = [vars sprintf('x%u,',i)];
	end
	vars = [vars sprintf('x%u]',nvars(eq(1)))];
		
	s = '[';
	for i=1:numel(eq)-1
		s = [s sprintf('poly(%s,%s,Dom::IntegerMod(%d)),',char(sym(eq(i),0)),vars,zp)];
	end
	s = [s sprintf('poly(%s,%s,Dom::IntegerMod(%d))]',char(sym(eq(end),0)),vars,zp)];
	
end