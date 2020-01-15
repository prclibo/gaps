classdef prob_pc_relpose_4pra_t2d < problem
    properties
        N = 4;
        q = sym('q%d%d', [3, 4]);
        qq = sym('qq%d%d', [3, 4]);
        u = sym('u%d', [3, 1]);
        s = sym('s');
        R;
    end
    methods
        function [eqs_sym, abbr_subs, unk_vars] = gen_eqs_sym(obj)
            obj.R = sym_quat2dcm([obj.s; obj.u]);
            
            eqs_sym = sym([]);
            abbr_subs = struct();
            for index = 1:obj.N
                [F, F_subs] = obj.Fijk(index);
                eqs_sym(end + 1) = det(F);
                abbr_subs = catstruct(abbr_subs, F_subs);
            end
            eqs_sym(end + 1) = obj.u.' * obj.u + obj.s .* obj.s - 1;
            unk_vars = [obj.u;].';
        end
        function [in_subs, out_subs] = gen_par_subs(obj)
            in_subs = struct();
            in_subs.q = obj.q;
            in_subs.qq = obj.qq;
            in_subs.s = obj.s;
            out_subs.u = obj.u;
        end
        function [in_zp, out_zp] = rand_var_zp(obj, p)
            [in_zp, out_zp] = vzp_rand_relpose_pinhole(4, p, {'q', 'qq', 's'}, {'u'});
        end
        function [in_rl, out_rl] = rand_par_rl(obj)
            [in_rl, out_rl] = vrl_rand_relpose_pinhole(4, {'q', 'qq', 's'}, {'u'});
        end
        function [F, abbr_subs] = Fijk(obj, i)
            j = mod(i, obj.N) + 1;
            k = mod(i + 1, obj.N) + 1;
            
            pij = sym(sprintf('p_%d_%d_%%d', i, j), [3, 1]);
            pik = sym(sprintf('p_%d_%d_%%d', i, k), [3, 1]);
            ppij = sym(sprintf('pp_%d_%d_%%d', i, j), [3, 1]);
            ppik = sym(sprintf('pp_%d_%d_%%d', i, k), [3, 1]);
            
            q = obj.q;
            qq = obj.qq;
            R = obj.R;
            F = [transpose(qq(:, j)) * R * pij, transpose(ppij) * R * q(:, j);
                transpose(qq(:, k)) * R * pik, transpose(ppik) * R * q(:, k)];
            pij_val = skew(q(:, i)) * q(:, j);
            pik_val = skew(q(:, i)) * q(:, k);
            ppij_val = skew(qq(:, i)) * qq(:, j);
            ppik_val = skew(qq(:, i)) * qq(:, k);
            
            c = num2cell([pij_val; pik_val; ppij_val; ppik_val]);
            f = cellstr(string([pij, pik, ppij, ppik]));
            abbr_subs = cell2struct(c, f);
            
        end
    end
end