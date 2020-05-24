function str = cg_extract_solutions(solv,template,opt)
str = [];

scheme = template.extraction_scheme;
prob = solv.prob;
unk_vars = prob.unk_vars;
assert(numel(unk_vars) == numel(scheme));
for k = 1:length(scheme)
    str = [str sprintf('sols(%d,:) = %s;\n',scheme{k}.target_idx,cg_parse_scheme(scheme{k}))];
end

for i = 1:numel(prob.unk_vars)
    str = [str sprintf('%s = sols(%d, :);\n', char(prob.unk_vars(i)), i)];
end

if length(scheme) < solv.n_vars
    str = [str sprintf('warning(''TODO: Extract remaining variables.'');\n')];
end

str = [str, sprintf('nsols = numel(%s);\n', unk_vars(1))];
par_names = fieldnames(prob.out_subs);
for i = 1:numel(par_names)
    par_name = par_names{i};
    subs_ = prob.out_subs.(par_name);
    if numel(subs_) > 1
        str = [str, sprintf('%s = cell(1, nsols);\n', par_name)];
        str = [str, sprintf('for isol = 1:nsols\n')];
        str = [str, sprintf('\t%s{isol} = nan(%d, %d);\n',...
            par_name, size(subs_, 1), size(subs_, 2))];
        for ivar = 1:numel(subs_)
            str = [str, sprintf('\t%s{isol}(%d) = %s(isol); \n', par_name, ivar, char(subs_(ivar)))];
        end
        str = [str, sprintf('end\n')];
    else
        str = [str, sprintf('%s = num2cell(%s);\n', par_name, par_name)];        
    end
end

% Unpack grouped known vars.

