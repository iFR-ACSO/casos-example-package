%% -----------------------------------------------------------------------
%
% Short description: This tutorial describes how to setup a simple convex
%                    optimization problem with DSOS decsision variable and 
%                    DSOS constraint. In a simple step by step
%                    tutorial we guide you through the process of setting
%                    up a problem in CaΣoS, solve the problem and evaluate
%                    the solution statistcs.
%
%
%  Problem: For a given Lypunov function, we want to check if 
%           {x \in R^n | V(x) < 1} is a stable level set. 
%            We define a convex SOS problem of the form:
%
%           s*(V-1) - \nabla V \cdot f(x) - epsilon(x) is SDSOS 
%
%           where s is a SDSOS multiplier coming from the S-procedure, V is 
%           a given Lyapunov function, f is the system dynamic and 
%           epsilon(x) is a small term to ensure strict negativity.
%
%           s is the only sought SDSOS decision variable. 
%           No cost is assigned and no parameter are needed.
%
%   License: see License file of repository
%
% -----------------------------------------------------------------------


%% Step 1) Define indeterminate 2 variables and system dynamics

% system states
x = casos.Indeterminates('x',2);

% system dynamics
f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)];

%% Step 2) We want to verify, if a given Lyapunov function is globally valid

% Lyapunov function candidate
V = 1.5*x(1)^2 - x(1)*x(2) + x(2)^2;

s = casos.PS.sym('s', monomials(x, 1), 'gram'); % SOS multiplier

% derivative w.r.t. time
Vdot = nabla(V,x)*f;

% enforce positivity
l = 1e-6*(x'*x);

%% Step 3) Setup up the problem struct

% Step 3.1: Define problem 

% polynomial decision variables
x_lin   = [];  % linear decision variables
x_sos   = [];  % sos decision variables
x_dsos  = [];  % dsos decision variables
x_sdsos = s;   % sdsos decision variables

% setup constraints
g_lin   = [];                % linear constraints
g_sos   = [];                % sos constraints
g_dsos  = [];                % dsos constraints
g_sdsos = s*(V-1)-Vdot - l;  % sdsos constraints


% cost function
f = [];

% parameter
p = [];

% Step 3.2: Setup the struct
sos = struct();

% in struct: linear constraint first, then sos, dsos, sdsos
sos.g = [g_lin;g_sos;g_dsos;g_sdsos];

% in struct: linear constraint first, then SOS, dsos, sdsos
sos.x = [x_lin;x_sos;x_dsos;x_sdsos ];

if ~isempty(f)
    sos.f = f;
else
    % f no cost, do not add to problem struct
end

% parameter
sos.p = p;

% Step 3.3: Provide the problem size i.e. size of cones
nx_sos   = length(x_sos);
nx_dsos  = length(x_dsos);
nx_sdsos = length(x_sdsos);
nx_lin   = length(x_lin);

ng_sos   = length(g_sos);
ng_dsos  = length(g_dsos);
ng_sdsos = length(g_sdsos);
ng_lin   = length(g_lin);

opts.Kx.sos    = nx_sos;
opts.Kx.dsos   = nx_dsos;
opts.Kx.sdsos  = nx_sdsos;

opts.Kx.lin   = nx_lin;
opts.Kc.sos   = ng_sos;
opts.Kc.dsos  = ng_dsos;
opts.Kc.sdsos = ng_sdsos;

% Step 3.4: Solver options

% if true, error returns infeasible, but we want to do our own check below
opts.error_on_fail = false; 
% turn off Newton polytop reduction
opts.newton_solver = [];    


%% Step 4) Generate a CaSoS solver instance

sdp_solver = 'mosek'; % 'mosek', 'scs', 'clarabel'

S = casos.sossol('S', ...        % name of solver
                 sdp_solver, ... % SDP solver
                 sos, ...        % problem structure
                 opts);          % options for solver


%% Step 5) Call the solver to solve the convex SDP

% call solver and store solution in solution struct 'sol'
sol = S();

% check the solution status with the unified solution status
if strcmp(S.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS')
    disp('Succesful!')
else
   disp('Unsuccesful!')
end

%% Step 6) Check the solver statistics and problem size of the conic problem
% check for solver statistics; output differs in the solver
switch sdp_solver
    case 'sedumi'
        fprintf('Sedumi needed %d iterations and it took %d seconds\n',[S.stats.iter,S.stats.cpusec])

    case 'mosek'
        fprintf('Mosek needed %d iterations and it took %d seconds\n',[S.stats.mosek_info.MSK_IINF_INTPNT_ITER,S.stats.mosek_info.MSK_DINF_INTPNT_TIME])
end

% check the problem size of the underlying conic problem
fprintf('The A matrix has %d rows and  %d columns\n',[ S.info.sdpsol.conic.size_a(1), S.info.sdpsol.conic.size_a(2)])


%% Step 7) Extracting the polynomial solution

% get the decision variable (here the s-multiplier)
sol.x
% get cost (zero because we did not assigned it)
sol.f
% the constraints evaluated at the solution
sol.g
% dual variables
sol.lam_g
sol.lam_x




