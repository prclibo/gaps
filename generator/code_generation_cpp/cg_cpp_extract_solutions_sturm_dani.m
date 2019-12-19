function str = cg_cpp_extract_solutions_sturm_dani(solv,template,opt)

%TODO: scheme = template.extraction_scheme;

basis = template.basis;
xx = create_vars(solv.n_vars);
ind = find_mon_indices(basis,[xx;1]);
assert(all(ind>0),'NYI: currently sturm_dani only supports cases when all variables (and unit) are available in basis.')
str = [];
for i = 1:solv.n_vars
    str = [str sprintf('%ssols(%d,i) = v(%d)/v(%d);\n',[opt.cg_indentation opt.cg_indentation],i-1,ind(i)-1,ind(end)-1)];
end


