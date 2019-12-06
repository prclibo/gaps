function [eqs ind] = geneq(eq,degree,multpolys)
% eqs = geneq(eqs,degree,multpolys) generates an expanded set of
% equations by multiplying with all monomials up to 'degree'.
% If 'multpolys' is supplied, these polynomials are used as a "basis"
% instead of the monomials x_1,..,x_n

eq = eqsize(eq);


if nargin<3
	p = nvars(eq(1));
	x(1) = multipol(1);
	for k=p:-1:1
		var = zeros(p,1);
		var(k) = 1;
		x(k+1) = multipol(1,var);
	end
else
	x = eqsize(multpolys);
end

mon = multipol(1);
for k=1:degree
	mon = unique([1; mon(:)]*x(:)')';
end
mon = sortmons(mon,'plex');
eqs = eq(:)*mon(:)';
eqs = eqs(:);

% disp(['Multiplying with ' num2str(numel(mon)) ' monomials:']);
% disp(char(mon));

if nargout>1
	[C M] = polynomials2matrix(eq,'plex');
	[i j v] = find(C.');
	v = 1:numel(v);
	C = sparse(j,i,v,size(C,1),size(C,2));
	ind_eq = C*M;
	ind_eqs = ind_eq(:)*mon(:)';
	C2 = polynomials2matrix(ind_eqs(:));
	[i j v] = find(C2);
	ind = [i(:)'; j(:)'; v(:)'];
end