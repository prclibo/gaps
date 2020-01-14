GAPS: A Generator for Automatic Polynomial Solvers
==================================================

GAPS is a tool to generate automatic polynomial solvers for a given multi-var polynomials system with varying coefficients. It is originally intended to construct solvers for _minimal problems in computer vision_.

GAPS wraps and improves the `autogen_v0_5` from Viktor Larsson.

Quick Start
-----------

To construct a polynomial solver, inherit the `problem` class to specify your polynomial system. Implement three functions in your inheritance.

`[eqs, unk_vars] = gen_eqs_sym(obj)` creates sym equation polynomials and sym unknown variables. 

`[kwn_zp, unk_zp] = rand_var_zp(obj, p)` generates random sample on Zp for variables in this problem. Field names in `kwn_zp` and `unk_zp` correspond to the known and unknown sym variables in the polynomials.

`[in_rl, out_rl] = rand_par_rl(obj)` generates random sample on real field for variables in this problem. You should instantiate this member function for your problem. Field names in `kwn_rl` and `unk_rl` correspond to the known and unknown sym variables in the polynomials.

See `problems/prob_*.m` for examples.

After construct your problem, call `generate_solver` to run GAPS. Usually you will want to setup some options per problem, see `problems/gen_*.m` for examples.

