function str = cg_cpp_declare_function_parameters(solv, declare)

prob = solv.prob;

fname = fieldnames(prob.in_subs);
in_decl = cell(1, numel(fname));
for i = 1:numel(fname)
    subs_ = prob.in_subs.(fname{i});
    if declare
        if numel(subs_) == 1
            in_decl{i} = sprintf('double %s', fname{i});
        else
            in_decl{i} = sprintf('Eigen::MatrixXd const& %s', fname{i});
        end
    else
        in_decl{i} = fname{i};
    end
end

fname = fieldnames(prob.out_subs);
out_decl = cell(1, numel(fname));
for i = 1:numel(fname)
    subs_ = prob.out_subs.(fname{i});
    if declare
        out_decl{i} = sprintf('std::vector<Eigen::MatrixXcd>* %s', fname{i});
    else
        out_decl{i} = sprintf('std::vector<double>* %s', fname{i});
    end
end

str = strjoin([in_decl, out_decl], ', ');