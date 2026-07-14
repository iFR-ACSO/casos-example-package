%% -----------------------------------------------------------------------
%
% Short description: This tutorial demonstrates how to setup a simple
%                    polynomial optimization problem using the moment
%                    hierarchy in CaSOS. We minimize a ratio between two
%                    polynomials, formulated as a moment program.
%
%  Problem: Measure program formulation:
%           minimize <f, mu> 
%           subject to mu is a nonnegative measure on R
%                      <mu, g> = 1 
%           
%
%   License: see License file of repository
%
% -----------------------------------------------------------------------

%% Step 1) Define indeterminate variable and problem data

% indeterminate variable
x = casos.Indeterminates('x');

% polynomial cost function
f = x^2 - x;   
g = x^2 + 2*x + 1;  

% truncation degree for moments
deg = 3;

%% Step 2) Setup measure program variables

% moments up to degree 'deg' as decision variables
mu = casos.PS.sym('mu', monomials(x, 0:deg), 'dual');

%% Step 3) Setup the problem struct

% Step 3.1: Define problem components

% measure decision variables (moments)
x_meas = mu.primalize;

% constraints:
% 1) mu is a probability measure: <g, mu> = 1
g_lin = mu.evaluate(g) - 1;  % linear constraint
  
% cost function: <f, mu>
f_cost = mu.evaluate(f);

% Step 3.2: Setup the struct
sos = struct();

% constraints: linear cone constraints
sos.g = g_lin;

% decision variables: measure variables
sos.x = x_meas;

if ~isempty(f_cost)
    sos.f = f_cost;
else
    % no cost, do not add to problem struct
end

% Step 3.3: Provide the problem size i.e. size of cones

nx_meas = length(x_meas);

ng_lin = length(g_lin);

% decision variable cones
opts.Kx.meas = nx_meas;   % measure decision variables

% constraint cones
opts.Kc.lin  = ng_lin;    % linear constraints

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
fprintf('Minimum value: %.12g\n', full(sol.f));
