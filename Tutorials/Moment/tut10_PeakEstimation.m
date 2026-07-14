%% -----------------------------------------------------------------------
%
% Short description: This tutorial demonstrates how to set up a measure
%                    program for peak estimation of polynomial dynamical
%                    systems using the moment hierarchy in CaSOS. 
%                    We estimate the minimum value of a polynomial function
%                    over trajectories of a dynamical system by formulating
%                    a measure program.
%
%  Problem: Maximum value of a state function p(x) attained along trajectories:
%
%           P* = sup_{t*, x_0} p(x(t*| x_0)
%                s.t. t* \in [0,T] 
%                     x_0 \in X_0
%                     \dot{x} = f(t, x(t)), \forall t \in [0, T]
% 
%  Measure program formulation for peak estimation:
%       minimize <p, mu_term>
%       subject to mu_init is a nonnegative measure on initial states x(0)
%                  mu_occu is a nonnegative measure on occupation measures (t,x(t))
%                  mu_term is a nonnegative measure on terminal states (t,x(t))
%                  Liouville equation: mu_term = mu_init + L_f*mu_occu
%                  supp(mu_init) in {x | gx(x) >= 0}
%                  supp(mu_occu) in {(t,x) | gt(t) >= 0}
%                  supp(mu_term) in {(x,t) | gt(t) >= 0}
%
%   Interpretation: The occupation measure mu_occu encodes the distribution
%                   of visited states over time. The objective <p, mu_term>
%                   minimizes the expected value of p over the terminal
%                   measure.
%
%   Measures: - mu_init: initial occupation measure (at t=0)
%             - mu_occu: occupation measure along trajectories (t>0)
%             - mu_term: terminal occupation measure (at final time)
%           
%   Note: The Liouville equation encodes the dynamics via adjoint 
%         infinitesimal generator L_f* (which is named liouville)
%
%
%   License: see License file of repository
%
% -----------------------------------------------------------------------

%% Step 1) Define indeterminate variable and problem data

% indeterminate variables
x   = casos.Indeterminates('x', 2, 1); 
t   = casos.Indeterminates('t', 1, 1); 


% dynamical system
f = [x(2); 
     -x(1)-x(2)+(1/3)*x(1)^3];

% polynomial function to minimize
p = x(2); 

% terminal time
T   = 5;    

% truncation degree for moments for mu_occu
deg = 6;    

% get degree for mu_init and mu_term
degv = 2*deg-f.maxdeg;
degv = degv + rem(degv,2);

% support of the measures
supp_mu_init  = 0.4^2 - (1.*x(1)-1.5)^2 - x(2)^2; 
supp_mu_term  = t*(T-1.*t); 
supp_mu_occu  = t*(T-1.*t); 

%% Step 2) Setup measure program variables

% initial measure variable 
mu_init = casos.PS.sym('mu_init',monomials(x,0:degv/2), 'dual');

% terminal measure variable 
mu_term = casos.PS.sym('mu_term',monomials([t;x],0:degv/2), 'dual');

% occupation measure variable 
mu_occu = casos.PS.sym('mu_occu',monomials([t;x],0:deg), 'dual');

%% Step 3) Setup the problem struct

% Step 3.1: Define problem components

% measure decision variables (moments)
x_meas = [mu_init.primalize; 
          mu_term.primalize; 
          mu_occu.primalize];

% constraints:
% 1) Liouville equation
% 2) <mu_init,1> = 1 (probability measure)
% 3) support constraint: supp(mu_init) in {x | supp_mu_init(x) >= 0}
% 4) support constraint: supp(mu_term) in {(t,x) | supp_mu_term(t) >= 0}
% 5) support constraint: supp(mu_occu) in {(t,x) | supp_mu_occu(t) >= 0}

v = monomials([x;t], 0:(2*deg-f.maxdeg));
liouville_eq = project(mu_term.primalize, v)-project(mu_init.primalize, v)-primalize(mu_occu.liouville(f,x,t));

g_lin = [liouville_eq;
         mu_init.evaluate(1)-1
         ];

g_meas = [primalize(mu_init.support(supp_mu_init));    
          primalize(mu_term.support(supp_mu_term));    
          primalize(mu_occu.support(supp_mu_occu))
          ];

% cost function: <mu_term, p>
f_cost = mu_term.evaluate(p);

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
opts.Kc.meas = ng_meas;   % nonnegative measure cone constraint

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

fprintf('Minimum value: %.12g\n', full(sol.f));

%% Step 8) Plot random trajectories

mpf = casos.toolboxes.to_multipoly(f);
mpx = casos.toolboxes.to_multipoly(1.*x);
mpg = casos.toolboxes.to_multipoly(supp_mu_init);
Nsample = 1e2;
xin = psample(-mpg, mpx, [1.5;0], Nsample);

for i=1:Nsample
    x0 = xin(:,i);
    [xtraj,xconv] = psim(mpf,mpx,x0,T);
    figure(1)
    hold on
    plot(xtraj{2}(:,1), xtraj{2}(:,2), 'k')
end
xlim([-0.5,2])
