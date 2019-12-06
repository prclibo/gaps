function d = det(A)
% MULTIPOL determinant function
if(size(A,1)~=size(A,2)); error('Determinant only defined for square matrices'); end;

% ps = perms(1:size(A,1));
% d = multipol();
% for k = 1 : size(ps, 1)
%     d = d + permsign(ps(k,:)) * prod(A(sub2ind(size(A),[1:size(A,2)],ps(k,:))));
% end

% Compute using recursion and ignore zero elements for speed.
d = rec_det(A);

	function d = rec_det(A)
		if numel(A)==1
			d = A;
		elseif all(size(A)==[2 2])
			d = A(1,1)*A(2,2)-A(1,2)*A(2,1);
		else
			d = 0;
			for c=1:size(A,2)
				if ~all(coeffs(A(1,c))==0) % zero polynomial
					d = d + (mod(c,2)*2-1) * A(1,c) * rec_det(A(2:end,[1:c-1 c+1:end]));
				end
			end
		end
	end

end