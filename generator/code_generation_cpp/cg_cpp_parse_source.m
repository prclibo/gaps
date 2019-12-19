function src_str = cg_cpp_parse_source(sources)

% Fix 17.10.2018: added .array()

% Generate expressions
if strcmp(sources.type,'action')
    %src_str = 'D.transpose()';
    src_str = 'D.transpose().array()';
elseif strcmp(sources.type,'basis')
    src_str = sprintf('V.row(%d).array()',sources.idx-1);
    %src_str = sprintf('V.row(%d)',sources.idx-1);
elseif strcmp(sources.type,'reducible')
    src_str = sprintf('(RR.row(%d)*V.transpose().matrix()).array()',sources.idx-1);
    %src_str = sprintf('(RR.row(%d)*V.transpose().matrix())',sources.idx-1);
    warning('Not tested');
elseif strcmp(sources.type,'target')
    %src_str = sprintf('sols.row(%d)',sources.idx-1);
    src_str = sprintf('sols.row(%d).array()',sources.idx-1);
else
    error('Unknown source type? Tell Viktor.');
end

