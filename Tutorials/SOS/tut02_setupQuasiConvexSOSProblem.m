%% -----------------------------------------------------------------------
%
% Short description: This tutorial describes how to setup a simple quasi-
%                    convex optimization problem. In a simple step by step
%                    tutorial we guide you through the process of setting
%                    up a problem in CaΣoS, solve the problem and evaluate
%                    the solution statistcs.
%
%
%  Problem: For a given Lypunov function, we want to check if 
%           {x \in R^n | V(x) < \gamma} is a stable level set. 
%            We define a convex SOS problem of the form:
%
%           s*(V-\gamma) - \nabla V \cdot f(x) - epsilon(x) is SOS 
%
%           where s is a SOS multiplier coming from the S-procedure, gamma
%           is the sought stable level set V is a given Lyapunov function, 
%           f is the system dynamic  and epsilon(x) is a small term to 
%           ensure strict negativity.
%
%           s is the only sought SOS decision varible and gamma the sought
%           stable level set. No parameter are needed.
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

% SOS multiplier
s = casos.PS.sym('s', monomials(x, 1), 'gram'); 

% level of stability
g = casos.PS.sym('g');

% derivative w.r.t. time
Vdot = nabla(V,x)*f;

% enforce positivity
l = 1e-6*(x'*x);

%% Step 3) Setup up the problem struct

% Step 3.1: Define problem 

% polynomial decision variables
x_lin   = []; % linear decision variables
x_sos   = s;  % sos decision variables

% setup constraints
g_lin = [];                % linear constraints
g_sos =  s*(V-g)-Vdot - l; % SOS constraints

% cost function
f = -g;

% parameter
p = [];

% Step 3.2: Setup the struct
sos = struct();

% in struct: linear constraint first, then SOS
sos.g = [g_lin;g_sos];

% in struct: linear constraint first, then SOS
sos.x = [x_lin;x_sos];

if ~isempty(f)
    sos.f = f;
else
    % f no cost, do not add to problem struct
end

% parameter
sos.p = p;

% Step 3.3: Provide the problem size i.e. size of cones

nx_sos = length(x_sos);
nx_lin = length(x_lin);

ng_sos = length(g_sos);
ng_lin = length(g_lin);

% we need to prepare the underlying convex SOS problem
sdp_solver = 'mosek'; % 'mosek', 'scs', 'clarabel'

% options
opts = struct('sossol',sdp_solver);

opts.Kx.sos = nx_sos;
opts.Kx.lin = nx_lin;

opts.Kc.lin = ng_lin;
opts.Kc.sos = ng_sos;

% Step 3.4: Solver options for underlying convex SOS problem

% if true, error returns infeasible, but we want to do our own check below
opts.sossol_options.error_on_fail = false; 
% turn off Newton polytop reduction
opts.sossol_options.newton_solver = [];    

% turn on verbosity of bisection
opts.verbose = true;

%% Step 4) Generate a CaSoS solver instance

S = casos.qcsossol('S1','bisection',sos,opts);



%% Step 5) Call the bisection interface and extract solution

% call solver and store solution in solution struct 'sol'
startTimeBisec = tic;
sol = S();
endTimeBisect = toc;

% check the solution status with the unified solution status
if strcmp(S.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS')
    fprintf('The maximum stable level set is: %f \n',full(sol.f))
else
   disp('Unsuccesful!')
end

%% Demostration how to extract the single iteratations

% The bisection is called n timer, that is, we solve n SDP problems. You 
% can look at each results individually as demonstrated next
 
fprintf('Mosek needed %d iterations and it took %d seconds for the last bisection.\n',[S.stats.iter{end}.mosek_info.MSK_IINF_INTPNT_ITER,S.stats.iter{end}.mosek_info.MSK_DINF_INTPNT_TIME])

% check the problem size of the underlying conic problem (same size for all
% problems)
fprintf('The A matrix has %d rows and  %d columns.\n',[S.stats.iter{end}.conic.size_A(1),S.stats.iter{end}.conic.size_A(2)])
fprintf('The SDP has  %d decision variables.\n',S.stats.iter{end}.conic.n_decVar)




