function str = cg_cpp_extract_solutions(solv,template,opt)

scheme = template.extraction_scheme;
prob = solv.prob;
unk_vars = prob.unk_vars;
assert(numel(unk_vars) == numel(scheme));

str = [];

for k = 1:length(scheme)
    str = [str sprintf('%ssols.row(%d) = %s;\n',opt.cg_indentation,scheme{k}.target_idx-1,cg_cpp_parse_scheme(scheme{k}))];
end

for i = 1:numel(prob.unk_vars)
    str = [str sprintf('Eigen::VectorXcd %s = sols.row(%d);\n', char(prob.unk_vars(i)), i - 1)];
end

if length(scheme) < solv.n_vars
    str = [str sprintf('#error TODO: Extract remaining variables.\n')];
end

str = [str, sprintf('int nsols = %s.size();\n', unk_vars(1))];
par_names = fieldnames(prob.out_subs);
for i = 1:numel(par_names)
    par_name = par_names{i};
    subs_ = prob.out_subs.(par_name);
    str = [str, sprintf('for (int isol = 0; isol < nsols; ++isol) {\n')];
    str = [str, sprintf('\tEigen::MatrixXcd _%s(%d, %d);\n',...
        par_name, size(subs_, 1), size(subs_, 2))];
    for ivar = 1:numel(subs_)
        str = [str, sprintf('\t_%s(%d) = %s(isol); \n', par_name, ivar - 1, char(subs_(ivar)))];
    end
    str = [str, sprintf('\t%s->push_back(_%s);\n', par_name, par_name)];
    str = [str, sprintf('}\n')];
end
