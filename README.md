GAPS: A Generator for Automatic Polynomial Solvers
==================================================

[![image](https://img.shields.io/github/license/prclibo/gaps)](https://github.com/prclibo/gaps/blob/master/LICENSE)
[![image](https://img.shields.io/badge/cs.CV-arXiv%3A2004.11765-B31B1B.svg)](https://arxiv.org/abs/2004.11765)

GAPS is a tool to generate automatic polynomial solvers for a given multi-var polynomials system with varying coefficients. It is originally intended to construct solvers for _minimal problems in computer vision_.

GAPS wraps and improves the `autogen_v0_5` from Viktor Larsson.

References
----------

Cite the following report ([arXiv link](https://arxiv.org/abs/2004.11765)) if you are using GAPS. This report also includes a short tutorial of the software.

```bibtex
@article{li2020gaps,
  title={GAPS: Generator for Automatic Polynomial Solvers},
  author={Li, Bo and Larsson, Viktor},
  journal={arXiv preprint arXiv:2004.11765},
  year={2020}
}
```

Dependencies
------------

* GAPS is developed under Matlab 2018a. Feel free to open a pull request if you meet and fix compatibility issue.

* GAPS calls [`Macaulay2`](http://www.math.uiuc.edu/Macaulay2/) (GPLv2) to compute Groebner basis and requires it to be installed beforehand.

* GAPS contains a modified version of [LundUniversityComputerVision/multipol](https://github.com/LundUniversityComputerVision/multipol). The modification allows Matlab symbols to be coefficients of a multipol polynomial. We thank @kalleastrom for issuing a BSD 3-Clause license for multipol.


Quick Start
-----------

For longer tutorial and technique introdution, please check the [arXiv report](https://arxiv.org/abs/2004.11765).

To construct a polynomial solver, inherit the `problem` (see [`generator/problem.m`]()) class to specify your polynomial system. Implement three functions in your inheritance.

`[in, out] = gen_arg_subs(obj)` creates two structs corresponding to input/output variables. Field names of the struct will be argument names used in the generated function. Field values are sym variables that will be used to denote polynomials.

`[eqs, abbr_subs] = gen_eqs_sym(obj)` creates sym equation polynomials. 

`[in_zp, out_zp] = rand_arg_zp(obj, p)` generates random sample on Zp for variables in this problem. Field names in `kwn_zp` and `unk_zp` correspond to the known and unknown sym variables in the polynomials.

`[in_rl, out_rl] = rand_arg_rl(obj)` generates random sample on real field for variables in this problem. You should instantiate this member function for your problem. Field names in `kwn_rl` and `unk_rl` correspond to the known and unknown sym variables in the polynomials.

Below is a simple example on solving the 5-pt problem:

```matlab
classdef prob_pc_relpose_5p_nulle_ne__simple < problem
    % A simple instance to construct a pinhole camera 5-pt relative pose
    % estimation problem solver.
    methods
        function [in, out] = gen_arg_subs(obj)
            % Each field in `in/out` will become an input/output
            % argument in the generated solver.
            %
            % `in.NE` is a 4x3x3 matrix made up of NEijk symbols. The
            % created solver will expect input argument NE to be 4x3x3
            % matrix and fill its element values to NEijk respectively.
            %
            % See your favorite 5-pt paper for details.
            %
            % Base vectors of the null space of the essential matrix.
            in.NE = sym('NE%d%d%d', [4, 3, 3]);
            % The weights
            out.w = sym('w%d', [3, 1]);
        end
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            [in, out] = gen_arg_subs(obj);
            NE = permute(in.NE, [2, 3, 1]);
            NE = reshape(NE, 9, 4);
            E = reshape(NE * [out.w; 1], 3, 3);
            eqs_sym = sym([]);
            
            % Construct polynomial system as symbolics
            eqs_sym(1) = det(E);
            Et = transpose(E);
            te = 2*(E*Et)*E - trace(E*Et)*E;
            eqs_sym(2:10) = te(:);
            
            % abbr_subs (Abbreviation substitution) is used to declare
            % intermediate variables when computing coefficients and their
            % expansions. This this very simple case we will not use it.
            abbr_subs = struct([]);
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            % For 5-pt problem, arbitrary value of NE always corresponds to
            % valid polynomial system. Therefore we just instantiate a Zp
            % case by random integer. However, usually this does not hold
            % for a general minimal problem. You need to consider the
            % geometry constraint during instantiation.
            %
            % in_zp/out_zp is expected to have same fields as of in_subs
            % and out_subs returned by `gen_arg_subs`. But the field values
            % are symbolic integers.
            in_zp.NE = sym(randi([1, p - 1], [4, 3, 3]));
            % out_zp can be omitted as it is not used right now.
            out_zp = struct();
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            error(['This is similar to `rand_arg_zp` with real field values. ',...
                'We are not using it here as it is for benchmarking']);
        end
    end
end
```

See `problems/prob_*.m` for more examples.

After construct your problem, call `generate_solver` to run GAPS. Usually you will want to setup some options per problem, see `problems/gen_*.m` for examples.

