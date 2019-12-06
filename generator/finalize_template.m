function [ template ] = finalize_template( solv, template, opt )
% [template]=finalize_template(solv,template)
% Computes template indices and other stuff needed for creating the solver.

% Figure out how to normalize the eigenvectors and extract solutions
[template.extraction_scheme,template.norm_scheme] = ...
            find_extraction_scheme(solv,template,opt);

mdegs = monomials(solv.coefficients.ind_eqs);
cidxs = coeffs(solv.coefficients.ind_eqs);
neqs = sum(cellfun(@(x) size(x, 2), template.monomials));
eqs_templ = struct('mdegs', cell(1, neqs), 'cidxs', cell(1, neqs));
offset = 0;
for k = 1:numel(mdegs)
    for i = 1:size(template.monomials{k}, 2)
        offset = offset + 1;
        eqs_templ(offset).mdegs = mdegs{k} + template.monomials{k}(:, i);
        eqs_templ(offset).cidxs = cidxs{k};
    end
end
[cc, mm] = coeff2mat({eqs_templ.cidxs}, {eqs_templ.mdegs});

basis = template.basis;
reducible = template.reducible;

if ~isempty(opt.saturate_mon)
    basis = basis * template.satmon;
    reducible = reducible * template.satmon;
end

% reorder monomials
ind_basis = find_mon_indices(mm, basis);
if any(ind_basis == -1)
    % basis monomials missing from expanded set (and not needed)
    % HACK: just add them here with zero coefficient.
    cc = [cc zeros(size(cc,1),nnz(ind_basis==-1))];
    mm = [mm;basis(ind_basis==-1)'];
    ind_basis = find_mon_indices(mm, basis);
end
ind_red = find_mon_indices(mm, reducible);
ind_excessive = setdiff(1:length(mm),[ind_basis ind_red]);
cc = cc(:,[ind_excessive ind_red ind_basis]);
mm = mm(:, [ind_excessive ind_red ind_basis]); % FIX (libo)

template.C_sz = size(cc);
template.C_ind = find(cc);
[ii,jj] = find(cc);
template.C_ind2 = [ii jj];
template.C_coeff = cc(template.C_ind);
template.mm = mm;

C0 = cc(:,1:size(cc,2)-length(ind_basis));
C1 = cc(:,size(cc,2)-length(ind_basis)+1:end);

% Indices for template parts C = [C0 C1]
template.C0_sz = size(C0);
template.C1_sz = size(C1);
template.C0_ind = find(C0);
template.C1_ind = find(C1);
[ii,jj] = find(C0);
template.C0_ind2 = [ii jj];
[ii,jj] = find(C1);
template.C1_ind2 = [ii jj];
template.C0_coeff = C0(template.C0_ind);
template.C1_coeff = C1(template.C1_ind);

% Compute the indices of the reducible in the action matrix
% matrix RR maps basis to [reducible; basis]
available_mons = [template.reducible template.basis];
template.AM_ind = find_mon_indices(available_mons, template.action*template.basis);


