function [ templ ] = finalize_template( solv, templ, opt )
% [template]=finalize_template(solv,template)
% Computes template indices and other stuff needed for creating the solver.

% Figure out how to normalize the eigenvectors and extract solutions
[templ.extraction_scheme,templ.norm_scheme] = ...
            find_extraction_scheme(solv,templ,opt);

mdegs = monomials(solv.coefficients.ind_eqs);
cinds = coeffs(solv.coefficients.ind_eqs);
neqs = sum(cellfun(@(x) size(x, 2), templ.monomials));
eqs_templ = struct('mdegs', cell(1, neqs), 'cinds', cell(1, neqs));
offset = 0;
for k = 1:numel(mdegs)
    for i = 1:size(templ.monomials{k}, 2)
        offset = offset + 1;
        eqs_templ(offset).mdegs = mdegs{k} + templ.monomials{k}(:, i);
        eqs_templ(offset).cinds = cinds{k};
    end
end
[cc, mm] = coeff2mat({eqs_templ.cinds}, {eqs_templ.mdegs});

basis = templ.basis;
reducible = templ.reducible;

if ~isempty(opt.saturate_mon)
    basis = basis * templ.satmon;
    reducible = reducible * templ.satmon;
end

% reorder monomials
ind_bas = find_mon_indices(mm, basis);
if any(ind_bas == -1)
    % basis monomials missing from expanded set (and not needed)
    % HACK: just add them here with zero coefficient.
    cc = [cc zeros(size(cc,1),nnz(ind_bas==-1))];
    mm = [mm;basis(ind_bas==-1)'];
    ind_bas = find_mon_indices(mm, basis);
end
ind_red = find_mon_indices(mm, reducible);
ind_exc = setdiff(1:length(mm),[ind_bas ind_red]);
cc = cc(:,[ind_exc ind_red ind_bas]);
mm = mm(:, [ind_exc ind_red ind_bas]); % FIX (libo)

if opt.find_upper_trianglar
    assert(~opt.sparse_template);
    c_lm = arrayfun(@(r) find(cc(r, :), 1, 'first'), 1:size(cc, 1));
    [c_utri, r_utri] = unique(c_lm);
    
    row_perm = [r_utri(:).', setdiff(1:size(cc, 1), r_utri)];
    col_perm = [c_utri(:).', setdiff(1:size(cc, 2), c_utri)];
    cc = cc(row_perm, col_perm);
    mm = mm(:, col_perm);
    len_utri = numel(r_utri);
else
    len_utri = 0;
end

assert(len_utri <= numel(ind_exc));

rows = [len_utri, size(cc, 1) - len_utri];
cols = [len_utri, size(cc, 2) - len_utri - numel(ind_bas), numel(ind_bas)];
cind_blk = mat2cell(cc, rows, cols);
templ.cind_blk = cind_blk;
templ.cind_blk_rows = rows;
templ.cind_blk_cols = cols;

templ.C_sz = size(cc);
templ.C_ind = find(cc);
[ii,jj] = find(cc);
templ.C_ind2 = [ii jj];
templ.C_coeff = cc(templ.C_ind);
templ.mm = mm;

C0 = cc(:,1:size(cc,2)-length(ind_bas));
C1 = cc(:,size(cc,2)-length(ind_bas)+1:end);

% Indices for template parts C = [C0 C1]
templ.C0_sz = size(C0);
templ.C1_sz = size(C1);
templ.C0_ind = find(C0);
templ.C1_ind = find(C1);
[ii,jj] = find(C0);
templ.C0_ind2 = [ii jj];
[ii,jj] = find(C1);
templ.C1_ind2 = [ii jj];
templ.C0_coeff = C0(templ.C0_ind);
templ.C1_coeff = C1(templ.C1_ind);

% Compute the indices of the reducible in the action matrix
% matrix RR maps basis to [reducible; basis]
available_mons = [templ.reducible templ.basis];
templ.AM_ind = find_mon_indices(available_mons, templ.action*templ.basis);


