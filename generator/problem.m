classdef problem
    properties
        eqs = sym([]);
        eqs_sym = sym([]);
        
        abbr_subs = struct();
        unk_vars = sym([]);
        kwn_vars = sym([]);
        in_subs = struct();
        out_subs = struct();
    end
    methods
        function obj = problem()
            [obj.eqs_sym, obj.abbr_subs, obj.unk_vars] = obj.gen_eqs_sym();
            obj.eqs = multipol(obj.eqs_sym, 'vars', obj.unk_vars);

            all_vars = symvar(obj.eqs_sym);
            obj.kwn_vars = setdiff(all_vars, obj.unk_vars);
            
            % Check if reserved symbols are used
            matches = regexp(sym2char(all_vars), '\<C\d+\>', 'once');
            if ~isempty(matches)
                error('Found reserved symbols \<C\d+\>. Avoid using this!');
            end
            obj = obj.gen_coeff_subs();
            
            [obj.in_subs, obj.out_subs] = obj.gen_par_subs();
        end
        function [in_subs, out_subs] = gen_par_subs(obj)
            in_subs = cell2struct(num2cell(obj.kwn_vars), sym2charcell(obj.kwn_vars));
            out_subs = cell2struct(num2cell(obj.unk_vars), sym2charcell(obj.unk_vars));
        end
        function obj = gen_coeff_subs(obj)
            coeff_cell = cell(1, numel(obj.eqs));
            names = cell(1, numel(obj.eqs));
            idx = 0;
            for i = 1:numel(obj.eqs)
                coeff_cell{i} = coeffs(obj.eqs(i));
                names{i} = strseq('C', idx + (1:numel(coeff_cell{i})));
                names{i} = names{i}(:).';
                obj.eqs(i) = multipol(sym(names{i}), monomials(obj.eqs(i)),...
                    'vars', vars(obj.eqs(i)));
                idx = idx + numel(coeff_cell{i});
            end
            coeff_eqs = num2cell(cat(2, coeff_cell{:}));
            coeff_names = cat(2, names{:});
            coeff_subs = cell2struct(coeff_eqs(:), coeff_names(:));
            obj.abbr_subs = catstruct(obj.abbr_subs, coeff_subs);
        end
        function [eq_zp, unk_zp] = rand_eq_zp(obj, p)
            [kwn_zp, unk_zp] = obj.rand_var_zp(p);
            eq_zp = subs_var(obj.eqs_sym, catstruct(kwn_zp, obj.abbr_subs),...
                'verbose', 'zp', p);
            eq_zp = multipol(eq_zp, 'vars', obj.unk_vars);
        end
        
    end
    methods (Abstract)
        [kwn_zp, unk_zp] = rand_var_zp(obj, p)
        [in_rl, out_rl] = rand_par_rl(obj)
        [eqs, unk_vars] = gen_eqs_sym(obj)
    end
    methods (Static)
        function var_subs = unpack_pars(par_subs, val_subs)
            function y = expd_cat_fields(x)
                y = struct2cell(x);
                y = cellfun(@(z) z(:), y, 'UniformOutput', false);
                y = cat(1, y{:});
            end
            pars = sym2charcell(expd_cat_fields(par_subs));
            vals = num2cell(expd_cat_fields(val_subs));
            var_subs = cell2struct(vals, pars);
        end
    end
end