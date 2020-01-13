function str = cg_cpp_extract_solutions(solv,template,opt)

scheme = template.extraction_scheme;
prob = solv.prob;
unk_vars = prob.unk_vars;
assert(numel(unk_vars) == numel(scheme));

str = [];

% str = [str sprintf('std::cerr << D << std::endl;\n')];
% str = [str sprintf('std::cerr << V.transpose() << std::endl;\n')];


for k = 1:length(scheme)
    str = [str sprintf('%ssols.row(%d) = %s;\n',opt.cg_indentation,scheme{k}.target_idx-1,cg_cpp_parse_scheme(scheme{k}))];
end

if length(scheme) < solv.n_vars
    str = [str sprintf('#error TODO: Extract remaining variables.\n')];
end

