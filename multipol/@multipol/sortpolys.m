function [eq ind] = sortpolys(eq)

ind = 1:numel(eq);
quicksort(1,numel(eq));
eq = eq(ind);

	function quicksort(a,b)
		if a<b
			p = floor((a+b)/2);
			ind([p b]) = ind([b p]);
			s = a;
			for i=a:b-1
				if eq(ind(i))>eq(ind(b))
					ind([s i]) = ind([i s]);
					s = s+1;
				end
			end
			ind([s b]) = ind([b s]);
			quicksort(a,s-1);
			quicksort(s+1,b);
		end
	end
end