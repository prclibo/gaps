function display(s)
% MULTIPOL/DISPLAY Command window display of a multivariate polynomial object;
%

if isequal(get(0,'FormatSpacing'),'compact')
	sp = '';
else
	sp = ' ';
end

if numel(s)==1
	disp(sp);
	disp([inputname(1),' = '])
    deg_tab = array2table(s.monomials);
    deg_tab.Properties.VariableNames = cellstr(strcat('C', string(1:size(s.monomials, 2))));
    deg_tab.Properties.RowNames = cellstr(string(s.vars));
    if isempty(s.vars)
        disp(s.coeffs);
    else
        disp(deg_tab);
    end
elseif numel(s)>50 || isempty(s)
	disp(sp);
	disp([inputname(1),' = '])
	disp(sp);
	siz = size(s);
	fprintf('   %u',siz(1));
	for i=2:length(siz)
		fprintf(' x %u',siz(i));
	end
	fprintf(' multipol object\n')
	disp(sp);
else
	disp(sp);
	disp([inputname(1),' = '])
	disp(sp);
	disp(char(s));
	disp(sp);
end