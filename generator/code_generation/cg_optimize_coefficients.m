function [ coeff_str ] = cg_optimize_coefficients( coeffs_eq )
% Does essentially the same thing as 
% https://awful.codeplex.com/SourceControl/latest#matlab/au_ccode.m
% but generates Matlab code instead.

eq_sym = sym(coeffs_eq(:));
%eq_sym = simplify(eq_sym);
eq_sym = feval(symengine, 'generate::optimize', eq_sym);
c = feval(symengine, 'generate::MATLAB', eq_sym);

% Convert x1 -> data(1) and so on
coeff_str = char(c);
% for k = nvars(coeffs_eq(1)):-1:1
%     coeff_str = strrep(coeff_str,sprintf('x%d',k),sprintf('data(%d)',k));
% end

coeff_str = strrep(coeff_str,'t1 = ','coeffs = ');
coeff_str = strrep(coeff_str,'t1(','coeffs(');

% Workaround for MATLAB 2017b (and probably others)
coeff_str = strrep(coeff_str,'t1 = [','coeffs = [');


% hack to convert \n to actual newlines;
coeff_str = sprintf(coeff_str);
end

