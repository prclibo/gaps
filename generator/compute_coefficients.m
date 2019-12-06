function [ coefficients ] = compute_coefficients( prob, solv, opt )

generic_coeffs = opt.generic_coefficients;


if generic_coeffs
    fprintf('Using generic coefficients.\n');

    eqs = solv.eqs_zp{1};
    [cc,mm] = polynomials2matrix(eqs);
    cc = cc';
    ind = find(cc);
    cc(ind) = 1:length(ind);    
    cc = cc';
    ind_eqs = cc*mm;
    coeff_eqs = create_vars(length(ind));
else
    eqs = prob.eqs;
    eqs_zp = solv.eqs_zp;
    
    coeff_eqs = [];
    ind_eqs = [];
    coeff_zps = cell(1, numel(eqs_zp));
    for k = 1:solv.n_eqs
%         [coeff_eq, mdegs] = coeff2mat(coeffs(eqs(k)), monomials(eqs(k)));
        coeff_eq = coeffs(eqs(k));
        mdegs = monomials(eqs(k));
        for i = 1:numel(eqs_zp)
            coeff_zp = collect_coeffs(eqs_zp{i}(k), mdegs);
            coeff_zps{i} = [coeff_zps{i}, coeff_zp];
        end

        ind = length(coeff_eqs)+(1:+length(coeff_eq));
        ind_eqs = [ind_eqs; multipol(ind, mdegs, 'vars', vars(eqs(k)))];
        coeff_eqs = [coeff_eqs coeff_eq];        
    end
%     coeff_eqs = remove_vars(coeff_eqs,1:solv.n_vars);
end

coefficients.coeff_eqs = coeff_eqs;
coefficients.coeff_zps = coeff_zps;
coefficients.ind_eqs = ind_eqs;
% coefficients = prune_duplicate(coefficients);
end

function coefficients = prune_duplicate(coefficients)
[c,ia,ic] = unique(coefficients.coeff_eqs);
if length(c) == coefficients.coeff_eqs
    return;
end

% update indices
for k = 1:length(coefficients.ind_eqs)
    [cc,mm] = polynomials2matrix(coefficients.ind_eqs(k));
    coefficients.ind_eqs(k) = m2p(ic(cc),mm);
end
coefficients.coeff_eqs = c;
end
