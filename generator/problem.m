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
            [obj.eqs_sym, obj.abbr_subs] = obj.gen_eqs_sym();
            [obj.in_subs, obj.out_subs] = obj.gen_arg_subs();
            function out = flatten_and_cat_fields(in)
                out = struct2cell(in);
                out = cellfun(@(x) x(:).', out, 'UniformOutput', false);
                out = horzcat(out{:});
            end
            obj.kwn_vars = flatten_and_cat_fields(obj.in_subs);
            obj.unk_vars = flatten_and_cat_fields(obj.out_subs);            
            obj.eqs = multipol(obj.eqs_sym, 'vars', obj.unk_vars);
            
            assert(numel(unique(obj.kwn_vars)) == numel(obj.kwn_vars));
            assert(numel(unique(obj.unk_vars)) == numel(obj.unk_vars));
            
            all_vars = [obj.kwn_vars, obj.unk_vars];
            % Check if reserved symbols are used
            matches = regexp(sym2char(all_vars), '\<C\d+\>', 'once');
            if ~isempty(matches)
                error('Found reserved symbols \<C\d+\>. Avoid using this!');
            end
        end
        function [in, out] = gen_arg_subs(obj)
            in = cell2struct(num2cell(obj.kwn_vars), sym2charcell(obj.kwn_vars));
            out = cell2struct(num2cell(obj.unk_vars), sym2charcell(obj.unk_vars));
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
        function eq_zp = rand_eq_zp(obj, p)
            [in_zp, ~] = obj.rand_var_zp(p);
            kwn_zp = obj.unpack_pars(obj.in_subs, in_zp);
            
            eq_zp = subs_var(obj.eqs_sym, catstruct(kwn_zp, obj.abbr_subs),...
                'verbose', 'zp', p);
            eq_zp = multipol(eq_zp, 'vars', obj.unk_vars);
        end
        
        %------------------------------------------------------------------
        function [eqs, unk_vars] = gen_eqs_sym(obj)
            % gen_eqs_sym Creates polynomials
            %   [eqs, unk_vars] = gen_eqs_sym(obj) creates sym equation
            %   polynomials and unknown sym variables. You should instantiate
            %   this member function for your problem.
            error(['rand_var_zp is not implemented yet!',...
                'You should implement this according to your problem']);
        end
        %------------------------------------------------------------------
        function [in_zp, out_zp] = rand_var_zp(obj, p)
            % rand_var_zp Randomize variables from Zp
            %   [kwn_zp, unk_zp] = rand_var_zp(obj, p) generates random
            %   sample on Zp for variables in this problem. You should instantiate
            %   this member function for your problem. Field names in
            %   kwn_zp and unk_zp correspond to the known and unknown sym
            %   variables in the polynomials.
            
            error(['rand_var_zp is not implemented yet!',...
                'You should implement this according to your problem']);
        end
        %------------------------------------------------------------------
        function [in_rl, out_rl] = rand_par_rl(obj)
            % rand_var_rl Randomize variables from real field
            %   [kwn_rl, unk_rl] = rand_var_rl(obj) generates random
            %   sample on real field for variables in this problem. You should instantiate
            %   this member function for your problem. Field names in
            %   kwn_rl and unk_rl correspond to the known and unknown sym
            %   variables in the polynomials.
            
            error(['rand_var_rl is not implemented yet!',...
                'You should implement this according to your problem']);
            
        end
    end
    methods (Abstract)
    end
    methods (Static)
        %------------------------------------------------------------------
        function var_subs = unpack_pars(par_subs, val_subs)
            % unpack_pars unpacks sym elements as struct fields and sets values
            %   var_subs = unpack_pars(par_subs, val_subs)
            %
            %   Example
            %   -------
            %   % par_subs stores sym variables as fields like:
            %   par_subs.a = [a1, a2, a3];
            %   par_subs.b = [b1, b2; b3, b4];
            %
            %   % val_subs stores values (numeric or sym) as fields like:
            %   val_subs.a = [1, 2, 3];
            %   val_subs.b = [4, 5, 6, 7];
            %
            %   % var_subs will store values as unpacked sym variable names:
            %   % var_subs.a1 = 1;
            %   % var_subs.a2 = 2;
            %   % var_subs.a3 = 3;
            %   % var_subs.b1 = 4;
            %   % var_subs.b2 = 5;
            %   % var_subs.b3 = 6;
            %   % var_subs.b4 = 7;
            %   var_subs = unpack_pars(par_subs, val_subs)
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