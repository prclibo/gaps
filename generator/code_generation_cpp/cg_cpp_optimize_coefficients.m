function coeff_str = cg_cpp_optimize_coefficients(coeffs_eq)
% Does essentially the same thing as 
% https://awful.codeplex.com/SourceControl/latest#matlab/au_ccode.m
% but generates Matlab code instead which is converted to C++.

eq_sym = sym(coeffs_eq);
%eq_sym = simplify(eq_sym);
eq_sym = feval(symengine, 'generate::optimize', eq_sym);
c = feval(symengine, 'generate::C', eq_sym);


coeff_str = char(c);

coeff_str = strrep(coeff_str,'t1[0]','coeffs');
coeff_str = regexprep(coeff_str,'t(\d+)','_t$1_');
coeff_str = regexprep(coeff_str,'_t(\d+)_ = ','double _t$1_ = ');


coeff_str = sprintf(coeff_str);

return;




% 
% 
% eq_sym = sym(coeffs_eq);
% %eq_sym = simplify(eq_sym);
% eq_sym = feval(symengine, 'generate::optimize', eq_sym);
% c = feval(symengine, 'generate::MATLAB', eq_sym);
% 
% % Convert x1 -> d[0] and so on.
% coeff_str = char(c);
% for k = nvars(coeffs_eq(1)):-1:1
%     coeff_str = strrep(coeff_str,sprintf('x%d',k),sprintf('d[%d]',k-1));
% end
% 
% coeff_str = regexprep(coeff_str,'t1 = \[(.*)\];','coeffs << $1;');
% if contains(coeff_str,'t1(')
%     % TODO: Not sure when this comes into play...
%     % coeff_str = strrep(coeff_str,'t1(','coeffs(');
%     error('Not implemented.');
% end
% 
% coeff_str = regexprep(coeff_str,'t(\d+) = ','double t$1 = ');
% coeff_str = strrep(coeff_str,',',',\n');
% coeff_str = regexprep(coeff_str,'(d\[\d+\])\^(\d+)','std::pow($1,$2)');
% 
% % Convert \n to newlines.
% coeff_str = [sprintf(coeff_str) newline];
% end
% 
