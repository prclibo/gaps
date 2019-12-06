function s = size(mp);
if(numel(mp)==1)
    s = size(mp.coeffs,2);
else
    s = size(mp);
end
