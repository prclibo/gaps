function str = cg_cpp_parse_scheme(scheme)
src_str = cell(1,length(scheme.sources));
for k = 1:length(scheme.sources)
    src_str{k} = cg_cpp_parse_source(scheme.sources{k});
end

if strcmp(scheme.method,'direct')
    str = src_str{1};
elseif strcmp(scheme.method,'sqrt')
    str = sprintf('%s.sqrt()',src_str{1});
elseif strcmp(scheme.method,'division')
    str = sprintf('%s / (',src_str{1});
    for kk = 2:length(src_str)
        str = [str sprintf('%s',src_str{kk})];
        if kk < length(src_str)
            str = [str sprintf('*')];
        end
    end
    str = [str ')'];
elseif strcmp(scheme.method,'eigen') || strcmp(scheme.method,'sqrt_eigen')
    str = ['('];
    for kk = 1:length(src_str)
        str = [str sprintf('%s',src_str{kk})];
        if kk < length(src_str)
            str = [str sprintf('*')];
        end
    end
    str = [str ')'];
else
    error('Unknown scheme method? Tell Viktor?');
end

