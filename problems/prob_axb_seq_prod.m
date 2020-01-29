classdef prob_axb_seq_prod < problem
    methods
        function [in, out] = gen_arg_subs(obj)
            in.a1 = sym('a1_%d', [4, 1]);
            in.b1 = sym('b1_%d', [4, 1]);
            in.a2 = sym('a2_%d', [4, 1]);
            in.b2 = sym('b2_%d', [4, 1]);
            out.x = sym('x%d', [4, 1]);
        end
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            [in, out] = obj.gen_arg_subs();
            circ_prod = ...
                hami2lmat(in.a1) * hami2lmat(out.x) * hami2lmat(in.b1) *...
                hami2lmat(in.a2) * hami2lmat(out.x) * in.b2;
            eqs_sym = circ_prod - [1; 0; 0; 0];
            eqs_sym(5) = out.x.' * out.x - 1;
            
            abbr_subs = struct([]);
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            in_zp.a1 = zp_rand_unit(4, p);
            in_zp.b1 = zp_rand_unit(4, p);
            in_zp.a2 = zp_rand_unit(4, p);
            in_zp.b2 = zp_rand_unit(4, p);
            out_zp = struct();
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            error('');
        end
    end
end