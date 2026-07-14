%% -----------------------------------------------------------------------
%
% Short description: This tutorial demonstrates how to setup a simple
%                    polynomial optimization problem using the moment
%                    hierarchy in CaSOS. We minimize a 2-dimensional
%                    polynomial.
%                    
%  Problem: Minimize f(x) = 4x1^2+x1x2-4x2^2-2.1x1^4+4x2^4+x1^6/3 
%
%           We formulate the problem as a measure program:
%           minimize <f, mu> subject to mu is a probability measure on R^2
%           
%           The function has six local minima, two of them being global minima.
%
%   License: see License file of repository
%
% -----------------------------------------------------------------------

%% Step 1) Define indeterminate variable and problem data

% indeterminate variable
x = casos.Indeterminates('x', 2, 1); 

% polynomial cost function
f = 4*x(1)^2+x(1)*x(2)-4*x(2)^2-2.1*x(1)^4+4*x(2)^4+x(1)^6/3;
    
% truncation degree for moments
deg = 3;

%% Step 2) Setup measure program variables

% moments up to degree 'deg' as decision variables
mu = casos.PS.sym('mu', monomials(x, 0:deg), 'dual');

%% Step 3) Setup the problem struct

% Step 3.1: Define problem components

% linear decision variables (none in this case)
x_lin = [];

% measure decision variables (moments)
x_meas = mu.primalize;

% constraints:
% 1) mu is a probability measure: <1, mu> = 1
g_lin = mu.evaluate(1) - 1;  % linear constraint: sum of moments at degree 0 = 1

% measure cone constraints are handled via opts.Kc.meas
g_meas = [];    % support constraint

% cost function: <f, mu>
f_cost = mu.evaluate(f);

% parameters (none)
p = [];

% Step 3.2: Setup the struct
sos = struct();

% constraints: linear first, then measure cone constraints
sos.g = [g_lin; g_meas];

% decision variables: linear first, then measure variables
sos.x = [x_lin; x_meas];

if ~isempty(f_cost)
    sos.f = f_cost;
else
    % no cost, do not add to problem struct
end

% parameter
sos.p = p;

% Step 3.3: Provide the problem size i.e. size of cones

nx_meas = length(x_meas);
nx_lin = length(x_lin);

ng_lin = length(g_lin);
ng_meas = length(g_meas);

% decision variable cones
opts.Kx.lin  = nx_lin;    % linear decision variables
opts.Kx.meas = nx_meas;   % measure decision variables

% constraint cones
opts.Kc.lin  = ng_lin;    % linear constraints
opts.Kc.meas = ng_meas;   % measure cone constraints

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
sol = S('lbx', -inf, 'ubx', +inf, 'lbg', 0, 'ubg', 0);

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
sol.x.disp
fprintf('Minimum value: %.12g\n', full(sol.f));
