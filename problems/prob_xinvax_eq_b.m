classdef prob_xinvax_eq_b < problem
    properties
    end
    methods
        function obj = prob_xinvax_eq_b(config)
            obj = obj@problem(config);
        end
        function [in, out] = gen_arg_subs(obj)
            in.a = sym('a%d', [4, 1]);
            in.b = sym('b%d', [4, 1]);
            out.x = sym('x%d', [4, 1]);
        end
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            [in, out] = obj.gen_arg_subs();
            
            sgn = diag([1, -1, -1, -1]);
            lhs = hami2lmat(out.x) * hami2lmat(in.a) * sgn * out.x;
            
            eqs_sym = lhs - in.b;
            eqs_sym(5) = out.x.' * out.x - 1;
            
            abbr_subs = struct([]);
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            in_zp.a = zp_rand_unit(4, p);
            out_zp.x = zp_rand_unit(4, p);
            sgn = mod(diag([1, -1, -1, -1]), p);

            in_zp.b = mod(hami2lmat(out_zp.x) * hami2lmat(in_zp.a) * sgn * out_zp.x, p);
            return;
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            error('a');
            deg = obj.config.deg;
            out_rl.x = normc(rand(4, 1));
            in_rl.a = normc(rand(4, deg));
            in_rl.b = normc(rand(4, deg - 1));
            circ_prod = 1;
            for i = 1:deg
                sgn = diag([1, [1, 1, 1] * in_rl.pwr(i)]);
                if i < deg
                    axb = hami2lmat(in_rl.a(:, i)) * hami2lmat(sgn * out_rl.x) * hami2lmat(in_rl.b(:, i));
                    circ_prod = circ_prod * axb;
                else
                    ax = hami2lmat(in_rl.a(:, i)) * hami2lmat(sgn * out_rl.x);
                    circ_prod = circ_prod * ax;
                end
            end
            circ_prod = circ_prod * [1; 0; 0; 0];
            in_rl.b(:, deg) = quatinv(circ_prod(:).');

        end
    end
end
