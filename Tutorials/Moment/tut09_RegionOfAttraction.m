%% -----------------------------------------------------------------------
%
% Short description: This tutorial demonstrates how to setup a measure
%                    program for region of attraction estimation using the
%                    moment hierarchy in CaSOS. We estimate the region of
%                    attraction of the Van-der-Pol oscillator by formulating
%                    the occupation measure framework as a measure program.
%
%  Problem: Measure program formulation for region of attraction:
%           minimize -<1, mu_init>
%           subject to mu_init is a nonnegative measure on initial states x(0)
%                      mu_occu is a nonnegative measure on occupation measures (t,x(t))
%                      mu_term is a nonnegative measure on terminal states (t,x(t))
%                      Liouville equation: 
%                           mu_term = mu_init + L_f*mu_occu
%                      mu_init + (Lebesgue measure on region) is a measure
%                      supp(mu_init) in {x | gx1(x) >= 0, gx2(x) >= 0}
%                      supp(mu_occu) in {(t,x) | gx1(x) >= 0, gx2(x) >= 0}
%                      supp(mu_term) in {(x,t) | gxt(x) >= 0, gt(t) >= 0}
%
%   Interpretation: This measure program computes the region of attraction
%                   for the Van-der-Pol oscillator:
%                   dot(x1) = -2x2
%                   dot(x2) = 0.8x1 + 10(x1^2 - 0.21)x2
%                   
%                   The optimal measure mu_init* approximates the Lebesgue
%                   measure restricted to initial conditions that remain
%                   within the bounded region [-1.2,1.2]^2 and enter the
%                   terminal set {x: ||x||^2 <= 0.01} within time t in [0,100].
%                   The objective -<1, mu_init> minimizes the negative total mass,
%                   equivalently maximizing <1, mu_init>, which represents the
%                   volume of initial conditions in the region of attraction.
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
x = casos.Indeterminates('x', 2, 1);            
t = casos.Indeterminates('t', 1, 1);            

% Van-der-Pol oscillator
f = [-2*x(2); 0.8*x(1)+10*(x(1)^2-0.21)*x(2)]; 

% bounds on the region of attraction
gx1 = -(1.*x(1)-1.2)*(1.*x(1)+1.2);    
gx2 = -(1.*x(2)-1.2)*(1.*x(2)+1.2);    

% terminal region
gxt = (1e-2)-x'*x;                       

% bounds on time [0,100]
gt  = -1.*t*(1.*t-1e2);              

% truncation degree for moments for mu_occu
deg = 5;    

% get degree for mu_init and mu_term
degv = 2*deg-f.maxdeg;
degv = degv + rem(degv,2);

% get lebesgue measure on [-1.2 1.2]^2
v = monomials(x, 0:degv);
v = casos.PS(v.to_vector);
lam = int(v, x, [-1.2 1.2; -1.2 1.2])'*v;

%% Step 2) Setup measure program variables

% moments up to degree 'deg' as decision variables
mu_occu  = casos.PS.sym('mu',  monomials([t;x], 0:deg), 'dual');

% moments up to degree 'degv' as decision variables
mu_init = casos.PS.sym('mu0', monomials(x, 0:degv/2), 'dual');
mu_term = casos.PS.sym('mut', monomials([t;x], 0:degv/2), 'dual');

%% Step 3) Setup the problem struct

% Step 3.1: Define problem components

% measure decision variables (moments)
x_meas = [mu_init.primalize; 
          mu_term.primalize; 
          mu_occu.primalize];

% constraints:
% 1) Liouville equation
% 2) lam-mu_init is a nonnegative measure
% 3) support constraint: supp(mu_init) in {x | gx1 >= 0}
% 4) support constraint: supp(mu_init) in {x | gx2 >= 0}
% 5) support constraint: supp(mu_occu) in {x | gx1 >= 0}
% 6) support constraint: supp(mu_occu) in {x | gx2 >= 0}
% 7) support constraint: supp(mu_term) in {x | gxt >= 0}
% 8) support constraint: supp(mu_term) in {t | gt >= 0}

v = monomials([x;t], 0:2*deg-f.maxdeg);
liouville_eq = project(mu_term.primalize, v)-project(mu_init.primalize, v)-primalize(mu_occu.liouville(f,x,t));

g_lin = liouville_eq; % linear constraint

g_meas = [lam-mu_init.primalize;
          primalize(mu_init.support(gx1));
          primalize(mu_init.support(gx2));
          primalize(mu_occu.support(gx1));
          primalize(mu_occu.support(gx2));
          primalize(mu_term.support(gxt));
          primalize(mu_term.support(gt))];  

% cost function: -<mu_init, 1>
f_cost = -mu_init.evaluate(1);

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

fprintf('Volume upperbound value: %.12g\n', full(sol.f));

% check if Multipoly Toolbox is available
if exist('mpvar', 'file')
    figure(1)
    hold on
    pcontour(casos.toolboxes.to_multipoly(-sol.lam_x(1)), 1, [-1.2 1.2 -1.2 1.2])
end
