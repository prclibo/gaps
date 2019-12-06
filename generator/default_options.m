
function opt = default_options(opt)
if nargin < 1 || isempty(opt)
    opt = struct();
end


default_opt = [];


% Path to Macaulay2
default_opt.M2_path = 'M2';


% Action monomial. Can either be specified as an index or a vector.
% i.e. actmon = 3  =>  x3 
%  and actmon = [2;1;0]  =>  x1^2*x2
% Empty => tries all variables as action and takes smallest template
default_opt.actmon = [];

% If actmon is empty, take all basis monomials as actions
% Other all variables are taken as actions.
default_opt.actions_from_basis = 0;

% If all templates found should be constructed or just the "best" one
default_opt.build_all_templates = 0;

% Automatically handle systems with complex coefficients
default_opt.transform_C_to_Zp = 1;

% The number of integer expansions to run (i.e. how many integer instances
% should we run use).
default_opt.integer_expansions = 1;

% Matrix of [d, integer_expansions] where each column contains the data
% which is to be used for generating the integer instances.
% Default: [] => Let problem file generate random data.
default_opt.custom_data = [];

% Linear solver for computing the action matrix
% The options are 'backslash' (default) | 'qr' | 'pinv' 
% For the QR solver rank is determined as sum(r/r(1) > 1e-8)
default_opt.linear_solver = 'backslash';

% Eigensolver
% The options are 'default' | 'eigs_only' | 'sturm' (C++ only)
% For sturm and eigs_only we do the fast structured back subst. to solve
% for the eigenvectors. Note this only returns the real solutions.
% 'sturm' uses danilevskys method to get the characteristic poly and then
% uses Hartleys implementation of the sturm sequences
default_opt.eigen_solver = 'default';

% Characteristic poly solver (only for eigen_solver == 'sturm')
% How to find the characteristic polynomial for sturm solvers
% Valid options are 'la_budde' | 'danilevsky' | 'danilevsky_piv' | 'krylov'
default_opt.charpoly_method = 'danilevsky_piv';

% Custom monomial basis for the quotient space. Note that this does not
% have to be an actual basis, any set of spanning monomials is fine.
% If no basis is specified, the normal set is used.
default_opt.custom_basis = [];

% If the elimination template should be sparse. Use for large problems.
default_opt.sparse_template = 0;

% Removes 1 from the monomial basis if zero is a solution.
default_opt.remove_zero_sol = 0;

% Ensure extra monomials are reduced.
default_opt.extra_reducible = [];

% Force variables (x1,x2,x3,...) to be among the reducibles
% Note that this ensures a solution can always be extracted
% Setting this to 1 is equivalent to extra_reducible = create_vars(nvars)
default_opt.force_vars_in_reducibles = 0;

% If the evaluation of the coefficients should be optimized using MuPAD
default_opt.optimize_coefficients = 0;

% If the template should be reduced using the Syzygy module.
default_opt.syzygy_reduction = 1;

% Remove unecessary equations after forming the template.
default_opt.equation_pruning = 0;

% If the solver should be written to file. (Useful for testing.)
default_opt.write_to_file = 1;

% Output level for M2 (momomommmomomo?)
default_opt.M2_gbTrace = 0; % 3

% Quit after generating the template. (Useful for testing.)
default_opt.stop_after_template = 0;

% Assigns one coefficient for each monomial in the equations
% This is also used if the problem does not provide equations with both
% data and unknowns (eqs_data)
default_opt.generic_coefficients = 0;

% Prime field to work in
default_opt.prime_field = 30097;

% Number of bits to store monomial exponents. Empty => default in M2
default_opt.M2_monomial_size = [];

% Order of variables
% (e.g. [4 3 2 1] to order the 4 variables in reverse order).
default_opt.variable_order = [];

% Monomial ordering in M2
% e.g. 'GRevLex', 'Lex', 'GLex'
default_opt.M2_monomial_ordering = 'GRevLex';

% Variable weights for the monomial ordering
% Can either be a string or a matrix of weights
% Note that the opt.variable_order is not taken into account!
default_opt.M2_weights = [];

% If symmetries should be removed.
default_opt.use_sym = 1;

% Check for variable aligned symmetries.
default_opt.find_sym = 1;

% The maximum p to check for find_sym=1.
default_opt.sym_max_p = 2;

% The maximum degree of actions to consider for symmetry.
default_opt.sym_max_action_deg = 2;


% These parameter specify the symmetries, sym_cc specify the weights
% and sym_pp the corresponding degree. If find_sym=1 these will be
% populated automatically. Their sizes should be
%   size(sym_cc) = (number of variables) x (number of symmetries)
%   size(sym_pp) = (number of symmetries) x 1
%
% e.g.
%    sym_cc = [1 0 0; 0 1 1]', sym_pp = [2;2]
% corresponds to two 2-fold symmetries, one in the x1 and one in (x2,x3)
default_opt.sym_cc = [];
default_opt.sym_pp = [];

% The remainder to use for the basis monomials. 
% This can either be specified by a single number, or one per symmetry.
default_opt.sym_rem = [];


% When generating the code we insert line breaks when printing
% coefficients. Keep default.
default_opt.max_str_len = 50;

% Use fast C++ code to read output from Macaulay2.
% Requires extract_monomials.cpp to be mexed.
default_opt.fast_monomial_extraction = 1;

% Try to remove redudant columns from the template
default_opt.remove_extra_columns = 1;

% if we should also try to remove some redundant rows
% (only happens if remove_extra_columns = 1)
default_opt.remove_extra_rows = 1;

default_opt.cg_language = 'matlab';
default_opt.cg_indentation = '    ';
default_opt.cg_compile_mex = 1;
default_opt.cg_eigen_dir = getenv('EIGEN_DIR');
% Maximum matrix size which should be allocated on stack, i.e.
default_opt.cg_max_stack_alloc_size = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimental stuff below
default_opt.saturate_mon = [];
default_opt.generalized_eigenvalue_solver = 0;

% merge into opt
fields = fieldnames(default_opt);
for k = 1:length(fields)
    if ~isfield(opt,fields{k}) || isempty(getfield(opt,fields{k}))
        opt = setfield(opt,fields{k},getfield(default_opt,fields{k}));
    end
end

opt = orderfields(opt);