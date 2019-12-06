function h = hash(p)
% MULTIPOL/HASH Compute a "hash" of the polynomial p for efficient set operations.

h = char(p,0,16);
