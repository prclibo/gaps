function m_vec = mons2vec(mons)
% m_vec = mons2vec(mons) takes all monomials in the multipol object mons
% and places them one by one in the multipol vector m_vec.
N = nterms(mons);
mm = monomials(mons);
for k = N:-1:1
    m_vec(k,1) = multipol(1, mm(:, k));
end