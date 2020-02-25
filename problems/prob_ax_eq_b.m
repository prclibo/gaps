classdef prob_ax_eq_b < problem
    methods
        function [in, out] = gen_arg_subs(obj)
            len = obj.config.len;
            in.a = sym('a%d%d', len);
            in.b = sym('b%d', [len, 1]);
            out.x = sym('x%d', [len, 1]);
        end
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            [in, out] = obj.gen_arg_subs();
            eqs_sym = in.a * out.x - in.b;
            abbr_subs = struct;
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            len = obj.config.len;
            in_zp.a = sym(randi(p, len) - 1);
            out_zp.x = sym(randi(p, [len, 1]) - 1);
            in_zp.b = sym(mod(in_zp.a * out_zp.x, p));
            
            RA = zp_rref(in_zp.a, p);
            assert(logical(RA(end) == 1));
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            len = obj.config.len;
            in_rl.a = rand(len);
            out_rl.x = rand([len, 1]);
            in_rl.b = in_rl.a * out_rl.x;
        end
    end
end