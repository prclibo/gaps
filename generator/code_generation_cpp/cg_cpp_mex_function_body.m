function str = cg_cpp_mex_function_body(solv, opt)

prob = solv.prob;

%     function [names, types] = extract_name_type(subs_)
%         
%         fname = fieldnames(subs_);
%         names = fieldnames(subs_);
%         mask = cellfun(@(x) numel(x) > 1, struct2cell(subs_));
%         types = {};
%         types(mask) = {'Eigen::MatrixXd'};
%         types(~mask) = {'double'};
%     end
% 
% [in_names, in_types] = extract_name_type(prob.in_subs);
% [out_names, out_types] = extract_name_type(prob.out_subs);
% 
% str = '';
% str = [str, strjoin(strcat(in_types(:), {' '}, in_names(:), ';'), '\n')];
% str = [str, strjoin(strcat(out_types(:), {' '}, out_names(:), ';'), '\n')];


in_name = fieldnames(prob.in_subs);
in_decl = cell(1, numel(in_name));
in_par = cell(1, numel(in_name));
for i = 1:numel(in_name)
    subs_ = prob.in_subs.(in_name{i});
    if numel(subs_) == 1
        in_decl{i} = sprintf('double %s = mxGetPr(prhs[%d])[0];\n', in_name{i}, i - 1);
    else
        in_decl{i} = sprintf('Eigen::Map<Eigen::MatrixXd const> %s(mxGetPr(prhs[%d]), %d, %d);\n',...
            in_name{i}, i - 1, size(subs_));
    end
    in_par{i} = sprintf('%s', in_name{i});
end

out_name = fieldnames(prob.out_subs);
out_decl = cell(1, numel(out_name));
out_par = cell(1, numel(out_name));
for i = 1:numel(out_name)
    subs_ = prob.out_subs.(out_name{i});
    if numel(subs_) == 1
        out_decl{i} = sprintf('std::vector<double> %s;\n', out_name{i});
    else
        out_decl{i} = sprintf('std::vector<Eigen::MatrixXcd> %s;\n', out_name{i});
    end
    out_par{i} = sprintf('&%s', out_name{i});
end

str = [opt.cg_indentation, strjoin([in_decl, out_decl], opt.cg_indentation)];

par_str = strjoin([in_par, out_par], ', ');
str = [str, sprintf('%ssolver_%s(%s);\n', opt.cg_indentation, solv.name, par_str)];


for i = 1:numel(out_name)
    str = [str, sprintf('%splhs[%d] = convertToMatlabCell(%s);\n', opt.cg_indentation, i - 1, out_name{i})];
end



