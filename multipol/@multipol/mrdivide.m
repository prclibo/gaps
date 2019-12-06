function [p3 a] = mrdivide(p1,p2)
% MULTIPOL/MRDIVIDE operator
% Division by scalar or polynomial division by set of multipols
% p3 = p1/p2

if isnumeric(p2)
	p1.coeffs = p1.coeffs/p2;
	p3 = p1;
	return
end

% Vad är det som händer? Hjälp
% else
%     monomials = fliplr(p1.monomials);
%     coeffs=fliplr(p1.coeffs);
%     max_coeff = max(coeffs);
%     p2.monomials = fliplr(p2.monomials);
%     p2.coeffs=fliplr(p2.coeffs);
%     for ii = 1:length(p1.coeffs)
%         for jj = 2:length(p2.coeffs)
%             mon_temp = monomials-(p2.monomials(:,jj)+monomials(:,ii))*ones(1,length(p1.coeffs));
%             col = find(sum(abs(mon_temp))==0);
%             if ~isempty(col)
%                 coeffs(col) = coeffs(col)-p2.coeffs(jj)*coeffs(ii)/p2.coeffs(1);
%                 coeffs(ii) = coeffs(ii)/p2.coeffs(1);
%                 if (abs(coeffs(col)/max_coeff))<(100*eps)
%                     coeffs(col) = 0;
%                 end
%             end
%         end
%     end
%     p3 = multipol(coeffs,monomials);
%     p3 = sort(p3);
%     p3 = squeeze(p3);
% end

if isnumeric(p1) && numel(p2)==1
	p3 = p1./p2;
else
	[p1 p2 zero] = eqsize(multipol(p1),multipol(p2),multipol);
	p3 = p1; % Output remainder
	
	for k=numel(p1):-1:1
		% Perform multivariate polynomial division
		p1(k) = sort(p1(k)); % grevlex order by default
		for i=1:numel(p2)
			p2(i) = sort(p2(i));
		end
		for i=numel(p2):-1:1
			a(k,i) = zero;
		end
		r = zero;
		id = p1(k);
		while id~=zero
			id = squeeze(id);
			LMid = id.monomials(:,1);
			for i=1:numel(p2) % Try every divisor in p2 in order until successful
				divisible = false;
				LMd = p2(i).monomials(:,1);
				if all(LMid-LMd>=0) % Divisible by p2(i)
					q = multipol(id.coeffs(1)/p2(i).coeffs(1),LMid-LMd);
					a(k,i) = a(k,i) + q;
					% Explicitly cancel leading term
					id = multipol(id.coeffs(2:end),id.monomials(:,2:end)) - multipol(p2(i).coeffs(2:end),p2(i).monomials(:,2:end))*q;
					divisible = true;
					break;
				end
			end
			if ~divisible
				% No division possible, move leading term to remainder and try again
				r = r + multipol(id.coeffs(1),id.monomials(:,1));
				id = multipol(id.coeffs(2:end),id.monomials(:,2:end));
			end
		end
		p3(k) = r;
	end
end

