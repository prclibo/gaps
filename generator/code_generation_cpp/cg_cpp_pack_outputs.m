function str = cg_cpp_pack_outputs(solv,template,opt)

prob = solv.prob;
unk_vars = prob.unk_vars;

str = [];
for i = 1:numel(prob.unk_vars)
    str = [str sprintf('Eigen::VectorXcd %s = sols.row(%d);\n', char(prob.unk_vars(i)), i - 1)];
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
