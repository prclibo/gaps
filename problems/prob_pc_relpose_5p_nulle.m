classdef prob_pc_relpose_5p_nulle < problem
    methods
        function obj = prob_pc_relpose_5p_nulle()
            E1 = sym('NE_1_%d_%d', [3 ,3]);
            E2 = sym('NE_2_%d_%d', [3 ,3]);
            E3 = sym('NE_3_%d_%d', [3 ,3]);
            E4 = sym('NE_4_%d_%d', [3 ,3]);
            syms w_1 w_2 w_3
            E = w_1*E1 + w_2*E2 + w_3*E3 + E4;
            
            eqs = sym(zeros(10, 1));
            eqs(1) = det(E);
            
            Et = transpose(E);
            te = 2*(E*Et)*E - trace(E*Et)*E;
            eqs(2:10) = te(:);
            
            unk_vars = [w_1, w_2, w_3];
            
            obj.eqs = eqs;
            obj.unk_vars = unk_vars;
        end
        function var_zp = rand_var_zp(obj, p)
            var_zp = vzp_rand_relpose_pinhole(5, p);
        end
    end
end
