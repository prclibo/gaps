classdef prob_simple_quadratic < problem
    methods
        function [in_subs, out_subs] = gen_arg_subs(obj)
            in_subs.c = sym('c%d%d', [2, 3]);
            out_subs.x = sym('x', [1, 2]);
        end
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            [in, out] = obj.gen_arg_subs();
            x = out.x(1); y = out.x(2); c = in.c;
            eqs_sym(1) = c(1, 1) * x^2 + c(1, 2) * x * y + c(1, 3);
            eqs_sym(2) = c(2, 1) * y^2 + c(2, 2) * x * y + c(2, 3);
            abbr_subs = struct;
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            in_zp.c = mod([     8     1   -38;    3     5     -57], p);
            out_zp.x = [2, 3];
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            in_rl.c = [     8     1   -38;    3     5     -57];
            out_rl.x = [2, 3];
        end
    end
    
end