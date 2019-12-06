function vars = get_m2_var(idxs)

vars = arrayfun(@(i) sym(sprintf('x%d', i)), idxs);