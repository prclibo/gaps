function [eqs ind] = geneq2(eq,degrees,multpolys)
% eqs = geneq(eqs,degrees,multpolys) generates an expanded set of
% equations by multiplying with all monomials with the property that the
% individual degree of variable x_i does not exceed degrees(i), and the
% total degree does not exceed max(degrees).
% If 'multpolys' is supplied, these polynomials are used as a "basis"
% instead of the monomials x_1,..,x_n (and the total degree may be greater
% than max(degrees).

eq = eqsize(eq);

if nargin<3
	x = [1; create_vars(nvars(eq(1)))];
else
	[eq,x] = eqsize(eq,multpolys);
end

if numel(degrees)~=nvars(eq(1))
	error('"degrees" must have as many elements as there are variables');
end


p = sum(x)^(max(degrees));
m = monomials(p);
m(:,any(bsxfun(@gt,m,degrees(:)))) = [];
p = multipol(ones(1,size(m,2)),m);
[~,mon] = polynomials2matrix(p,'plex');
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
