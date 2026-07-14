%% -----------------------------------------------------------------------
%
% Short description: This tutorial demonstrates how to setup a simple
%                    polynomial optimization problem using the moment
%                    hierarchy in CaSOS. We minimize a ratio between two
%                    polynomials, formulated as a moment program.
%
%  Problem: Measure program formulation:
%           minimize -<1, mu>
%           subject to mu is a nonnegative measure on R^2
%                      (Lebesgue measure on [-1,1]^2 - mu) is a nonnegative measure
%                      supp(mu) in K = {x | g(x) >= 0}
%
%  Interpretation: The optimal measure mu* approximates the Lebesgue 
%                   measure restricted to the set {x | g(x) >= 0} and [-1,1]^2,
%                   and <1, mu*> approximates its volume.
%
%
%   License: see License file of repository
%
% -----------------------------------------------------------------------

%% Step 1) Define indeterminate variable and problem data

% indeterminate variable
x = casos.Indeterminates('x', 2, 1);

% polynomial defining the basic semialgebraic set
g = x(1)*(x(1)^2 + x(2)^2)-(x(1)^4+x(1)^2*x(2)^2+x(2)^4);

% truncation degree for moments
deg = 5;

% get lebesgue measure moments on [-1 1]^2 up to deg
v = monomials(x, 0:2*deg);
v = casos.PS(v.to_vector);
lam = int(v, x, [-1 1; -1 1])'*v;

%% Step 2) Setup measure program variables

% moments up to degree 'deg' as decision variables
mu = casos.PS.sym('mu', monomials(x, 0:deg), 'dual');

%% Step 3) Setup the problem struct

% Step 3.1: Define problem components

% measure decision variables (moments)
x_meas = mu.primalize;

% constraints:
% 1) lam - mu is a nonnegative measure
% 2) support constraint: supp(mu) in K 
% measure cone constraints are handled via opts.Kc.meas
g_meas = [lam-mu.primalize;
          primalize(mu.support(g))];  
  
% cost function: -<f, mu>
f_cost = -mu.evaluate(1);

% Step 3.2: Setup the struct
sos = struct();

% constraints: nonnegative measure cone constraints
sos.g = g_meas;

% decision variables: measure variables
sos.x = x_meas;

% cost function to struct
sos.f = f_cost;

% Step 3.3: Provide the problem size i.e. size of cones
nx_meas = length(x_meas);
ng_meas = length(g_meas);

% decision variable cones
opts.Kx.meas = nx_meas;   % measure decision variables

% constraint cones
opts.Kc.meas  = ng_meas;    % linear constraints

% Step 3.4: Solver options

% if true, error returns infeasible
opts.error_on_fail = false;

%% Step 4) Generate a CaSoS solver instance

sdp_solver = 'mosek'; % 'mosek', 'scs', 'clarabel', 'sedumi'

fprintf('Generating solver instance using %s...\n', sdp_solver);
S = casos.sossol('S', ...        % name of solver
                 sdp_solver, ... % SDP solver
                 sos, ...        % problem structure
                 opts);          % options for solver

%% Step 5) Call the solver to solve the convex SDP

% call solver and store solution in solution struct 'sol'
sol = S('lbg', 0, 'ubg', 0);

% check the solution status with the unified solution status
if strcmp(S.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
    fprintf('Succesful!\n');
else
    fprintf('Unsuccesful!\n');
end

%% Step 6) Check the solver statistics and problem size of the conic problem
% check for solver statistics; output differs in the solver
switch sdp_solver
    case 'sedumi'
        fprintf('Sedumi needed %d iterations and it took %.5f seconds\n',[S.stats.iter,S.stats.cpusec])

    case 'mosek'
        fprintf('Mosek needed %d iterations and it took %.5f seconds\n',[S.stats.mosek_info.MSK_IINF_INTPNT_ITER,S.stats.mosek_info.MSK_DINF_INTPNT_TIME])
end

% check the problem size of the underlying conic problem
fprintf('The A matrix has %d rows and  %d columns\n',[S.stats.conic.size_A(1),S.stats.conic.size_A(2)])
fprintf('The SDP has  %d decision variables\n',S.stats.conic.n_decVar)

%% Step 7) Extracting the polynomial solution and display minimum value

fprintf('Moments polynomial:\n')
sol.x.remove_coeffs(1e-6).disp
fprintf('Volume upperbound value: %.12g\n', -full(sol.f));
