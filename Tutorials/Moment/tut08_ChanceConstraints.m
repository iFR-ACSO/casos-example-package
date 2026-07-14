%% -----------------------------------------------------------------------
%
% Short description: This tutorial demonstrates how to setup a simple
%                    polynomial optimization problem using the moment
%                    hierarchy in CaSOS. We compute the maximum probability 
%                    that a random variable w (uniform on [-1,1]) satisfies a
%                    polynomial inequality p(x,w) >= 0 for a deterministic 
%                    decision variable x in a basic semialgebraic set K, 
%                    formulated as a measure program.
%
%  Problem: Measure program formulation:
%           minimize -<1, mu> 
%           subject to mu is a nonnegative measure on (x,w)-space
%                      mu_x is a probability measure on x
%                      mu_x has support on K = {x | g(x) >= 0}
%                      mu has support on {(x,w) | p(x,w) >= 0}
%                      product measure between mu_x and the Lebesgue measure 
%                      on w is a probability minus mu must be a nonnegative measure
%
%  Interpretation: The optimal measure mu* approximates the joint 
%                   distribution of (x,w) restricted to p(x,w) >= 0,
%                   where w ~ Uniform[-1,1] and x is a deterministic
%                   parameter. The objective -<1, mu> minimizes the
%                   negative total mass, equivalently maximizing <1, mu>,
%                   which represents the probability P[p(x,w) >= 0].
%
%
%   License: see License file of repository
%
% -----------------------------------------------------------------------

%% Step 1) Define indeterminate variable and problem data

% indeterminate variable
x = casos.Indeterminates('x', 2, 1);  

% region of interest ({(x,w) | p(x,w)>=0})
p = 0.5*1.*x(2)*(1.*x(2)^2+(1.*x(1)-0.5)^2)-(1.*x(2)^4+1.*x(2)^2*(1.*x(1)-0.5)^2+(1.*x(1)-0.5)^4);            

% constraint on x
g = (1-1.*x(1))*(1.*x(1)+1);

% truncation degree for moments
deg = 5;  

% get moments for a uniform distribution over [-1, 1]
v = monomials(x(2), 0:2*deg);
v = casos.PS(v.to_vector);
mu_w = 0.5*int(v, x(2), [-1; 1])'*v;

%% Step 2) Setup measure program variables

% moments up to degree 'deg' as decision variables
mu   = casos.PS.sym('mu', monomials(x,0:deg), 'dual');
mu_x = casos.PS.sym('mu_x', monomials(x(1),0:deg), 'dual');

%% Step 3) Setup the problem struct

% Step 3.1: Define problem components

% measure decision variables (moments)
x_meas = [mu.primalize; 
          mu_x.primalize];

% constraints:
% 1) mu is a probability measure: <1, mu_x> = 1
% 2) prod(mu_x,mu_w) >= mu
% 3) support constraint: supp(mu_x) in {x | g(x) >= 0}
% 4) support constraint: supp(mu) in { (x,w) | p(x,w) >= 0}
g_lin = mu_x.evaluate(1) - 1;  % linear constraint

g_meas = [mu_x.primalize+mu_w-1-mu.primalize;
          primalize(mu_x.support(g));
          primalize(mu.support(p))];

% cost function: -<mu,1>
f_cost = -mu.evaluate(1);

% Step 3.2: Setup the struct
sos = struct();

% constraints: linear cone constraints
sos.g = [g_lin; g_meas];

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
ng_meas = length(g_meas);

% decision variable cones
opts.Kx.meas = nx_meas;   % measure decision variables

% constraint cones
opts.Kc.lin  = ng_lin;    % linear constraints
opts.Kc.meas = ng_meas;   % nonnegative measure cone constraints

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

fprintf('Probability upperbound value: %.12g\n', full(sol.f));
