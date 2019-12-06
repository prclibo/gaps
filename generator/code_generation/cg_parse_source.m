function src_str = cg_parse_source(sources)
% Generate expressions
if strcmp(sources.type,'action')
    src_str = 'diag(D).''';
elseif strcmp(sources.type,'basis')
    src_str = sprintf('V(%d,:)',sources.idx);
elseif strcmp(sources.type,'reducible');
    src_str = sprintf('RR(%d,:)*V',sources.idx);
elseif strcmp(sources.type,'target');
    src_str = sprintf('sols(%d,:)',sources.idx);
else
    error('Unknown source type? Tell Viktor.');
end