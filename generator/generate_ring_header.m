function [ str ] = generate_ring_header( nv, opt )

varstr = '';

order = 1:nv;
if ~isempty(opt.variable_order)
    order = order(opt.variable_order);
end

for i=order
    varstr = [varstr sprintf('x%u,',i)]; %#ok
end


str = sprintf('KK = ZZ / %d;\n',opt.prime_field);
extra = '';
if ~isempty(opt.M2_monomial_size)
    extra = sprintf(',MonomialSize=>%d',opt.M2_monomial_size);
end

if isfield(opt,'M2_monomial_ordering')
    extra = [extra sprintf(',MonomialOrder=>%s',opt.M2_monomial_ordering)];
end

if ~isempty(opt.M2_weights)
    w = opt.M2_weights;
    
    if isa(opt.M2_weights,'double')
        if size(w,2) == 1
            w = w';
        end
        extra = [extra ',MonomialOrder=>{'];
        for k = 1:size(w,1)
            extra = [extra 'Weights=>{' sprintf('%d,',w(k,:))];
            extra = [extra(1:end-1) '},'];
        end
        extra = [extra(1:end-1) '}'];
    else
        extra = [extra sprintf(',MonomialOrder=>{Weights=>{%s}}',w)];
    end
end

str = [str sprintf('R = KK[%s%s];\n',varstr(1:end-1),extra)];


end

