%% -----------------------------------------------------------------------
%
% Example 1: Aerostability of a cubesat with feathered geometry
%
% Short description: Verify almost global stability of a cubesat subject
%                    to aerodynamic torque in free molecular flow.
%
% For details see:
%
% [1] F. Geyer, F. Tuttas, W. Fichter, and T. Cunis,
% "Sum-of-Squares Stability Verification on Manifolds with Applications in
% Spacecraft Attitude Control," European Journal of Control, p. 101587,
% Jul. 2026, doi: 10.1016/j.ejcon.2026.101587.
%
% License: see License file of repository
%
% -----------------------------------------------------------------------

%% User settings

plot_flag = true;

%% Step 1) Define indeterminate variables and system dynamics

% State variables in original coordinates
xbar       = casos.PD('xb', 6);
xbar_star1 = [1; zeros(5,1)];     % stable equilibrium in original coordinates
xbar_star2 = [-1; zeros(5,1)];    % unstable equilibrium in original coordinates

% xbar(1:3): 2-axis attitude unit vector (negative wind direction)
% xbar(4:6): rotational rates expressed in body coordinates

% Aerodynamic torque (Sentman's method, polynomial fit) and rate damping
% are already contained in the generated dynamics.
fbar = example1_dynamics(xbar);

% Manifold constraint in original coordinates
hbar = xbar(1:3)'*xbar(1:3) - 1;


%% Step 2) Shift and scale system dynamics and define SOS quantities

% Coordinate transformation:
%
%   x = S*(xbar - xbar_star1)
%
S = diag([ones(3,1); 1/0.05.*ones(3,1)]);

% Transformed state
x = casos.PD('x', 6);

% Inverse map:
%
%   xbar = S\x + xbar_star1
%
% used to substitute the transformed state into the original expressions.
xbar_of_x = S\x + xbar_star1;

% Transform unstable equilibrium
xstar2 = S*(xbar_star2 - xbar_star1);

% Transform manifold constraint
h = subs(hbar, xbar, xbar_of_x);

% Transform dynamics:
%
%   xdot = S*fbar(inv(S)*x + xbar_star1)
%
f = S*subs(fbar, xbar, xbar_of_x);


% Unknown polynomial decision variables
p = casos.PS.sym('p', monomials(x, 0:6));         % equality constraint
V = casos.PS.sym('v', monomials(x, 1:2), 'gram'); % Lyapunov candidate


% Enforce positive definite V
eps1 = 1e-5;
l1 = eps1*(x'*x);

% Supply rate defined via equilibria
eps2 = 1e-5;
l2 = eps2*(x'*x)*((x-xstar2)'*(x-xstar2));


%% Step 3) Setup the problem struct

% Step 3.1: Define problem

% Polynomial decision variables
x_lin = p;     % linear polynomial decision variables
x_sos = V;     % SOS polynomial decision variables

% Setup constraints
g_lin = [];

% SOS relaxation of Theorem 1
g_sos = [
    V-l1;
    -nabla(V,x)*f + p'*h - l2
];

% Cost function
f_obj = 0;

% Parameters
par = [];


% Step 3.2: Setup the struct

sos = struct();

% In struct: linear constraints first, then SOS constraints
sos.g = [g_lin; g_sos];

% In struct: linear decision variables first, then SOS decision variables
sos.x = [x_lin; x_sos];

% Cost
if ~isempty(f_obj)
    sos.f = f_obj;
else
    % If no cost is assigned, do not add it to the problem struct.
end

% Parameters
sos.p = par;


% Step 3.3: Provide the problem size, i.e. size of cones

nx_sos = numel(x_sos);
nx_lin = numel(x_lin);

ng_sos = size(g_sos,1);
ng_lin = size(g_lin,1);

opts.Kx.sos = nx_sos;
opts.Kx.lin = nx_lin;

opts.Kc.lin = ng_lin;
opts.Kc.sos = ng_sos;


% Step 3.4: Solver options

% If true, an error is returned if the solver fails.
% Disable this so that the solver status can be checked below.
opts.error_on_fail = false;

% Turn off Newton polytope reduction / simplification
opts.newton_solver = [];


%% Step 4) Generate a CaSoS solver instance

sdp_solver = 'mosek'; % 'mosek', 'scs', 'clarabel'

Sos = casos.sossol( ...
    'S', ...          % name of solver
    sdp_solver, ...   % SDP solver
    sos, ...          % problem structure
    opts);            % solver options


%% Step 5) Call the solver to solve the convex SDP

% Call solver and store solution in solution struct 'sol'
sol = Sos();

% Check the solution status with the unified solution status
if strcmp(Sos.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
    disp('Successful!')
else
    disp('Unsuccessful!')
end


%% Step 6) Check solver statistics and problem size of the conic problem

% Check solver statistics; output differs depending on solver
switch sdp_solver

    case 'sedumi'
        fprintf( ...
            'Sedumi needed %d iterations and it took %d seconds\n', ...
            [Sos.stats.iter, Sos.stats.cpusec]);

    case 'mosek'
        fprintf( ...
            'Mosek needed %d iterations and it took %d seconds\n', ...
            [Sos.stats.mosek_info.MSK_IINF_INTPNT_ITER, ...
             Sos.stats.mosek_info.MSK_DINF_INTPNT_TIME]);

end


%% Step 7) Extract the polynomial solution

% Retrieve Lyapunov function solution and remove small coefficients
% that are numerically zero.
Vsol = remove_coeffs( ...
    sol.x(2), ...
    1e-10*max(sol.x(2).poly2basis));


%% Step 8) Create plot

if ~plot_flag
    return
end

% Lyapunov function in original coordinates
Vbar_sol = subs(Vsol, x, S*(xbar - xbar_star1));


% Initial condition used in paper
x0 = [
     0.9969
    -0.0792
    -0.0006
    -0.0359
     0.0304
    -0.0105
];


% Simulate the true, untransformed dynamics
[t,y] = ode45( ...
    @(~,xx) example1_dynamics(xx), ...
    linspace(0,800,501), ...
    x0);

t = t.';
y = y.'; % states in rows


% Evaluate Lyapunov function.
%
% to_function takes one scalar input per indeterminate, so pass the state
% trajectory row by row.
Vbar_solFun = Vbar_sol.to_function();

yrows = num2cell(y,2);

Vbar_val = full(Vbar_solFun(yrows{:}));


% Plot states and normalized Lyapunov function
figure('Name', 'States and Lyapunov function (original coordinates)');
clf;
hold on;

plot( ...
    t, ...
    y(1:3,:), ...
    LineStyle='--', ...
    LineWidth=.8);

plot( ...
    t, ...
    rad2deg(y(4:6,:)), ...
    LineWidth=.8);

plot( ...
    t, ...
    Vbar_val./max(Vbar_val), ...
    'k', ...
    LineWidth=1);

xlabel('$t \,[\mathrm{s}]$', 'Interpreter', 'latex');

legend( ...
    '$\bar{x}_1$', ...
    '$\bar{x}_2$', ...
    '$\bar{x}_3$', ...
    '$\bar{x}_4 \,[^{\circ}/s]$', ...
    '$\bar{x}_5 \,[^{\circ}/s]$', ...
    '$\bar{x}_6 \,[^{\circ}/s]$', ...
    '$\bar{V}_{\mathrm{s}}(\bar{x})$', ...
    'Interpreter', 'latex', ...
    'Location', 'SouthEast');

set( ...
    gca, ...
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 12);

hold off;


function f = example1_dynamics(in1)
%EXAMPLE1_DYNAMICS
%    F = EXAMPLE1_DYNAMICS(IN1)
%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    04-Aug-2026 09:28:10
%Aerostability of a feathered cubesat (Example 1).
%
%States:
%x(1:3)  2-axis attitude unit vector (negative wind direction)
%x(4:6)  body rates [rad/s]
%
%GENERATED FILE -- DO NOT EDIT.
%Regenerate with scripts/export_standalone.m in the paper repository
%https://github.com/Fabian-Geyer/ECC26-SOS-Verify-Manifold-Spacecraft
%
%Accepts either a numeric state vector or a casos.PD state vector.

xb1 = in1(1,:);
xb2 = in1(2,:);
xb3 = in1(3,:);
xb4 = in1(4,:);
xb5 = in1(5,:);
xb6 = in1(6,:);

f = [
    xb2.*xb6.*1.0 - xb3.*xb5.*1.0;
    xb1.*xb6.*-1.0 + xb3.*xb4.*1.0;
    xb1.*xb5.*1.0 - xb2.*xb4.*1.0;
    xb4.*-1.6e-2;
    xb3.*-1.131089808877158e-4 ...
        - xb5.*1.6e-2 ...
        + xb4.*xb6.*2.653061224489796e-1 ...
        - xb2.^2.*xb3.*4.036999974620292e-4 ...
        + xb2.^4.*xb3.*1.946161649467695e-4 ...
        - xb3.^3.*3.993936643812048e-4 ...
        + xb3.^5.*1.946161649467695e-4;
    xb2.*1.131089808877158e-4 ...
        - xb6.*1.6e-2 ...
        - xb4.*xb5.*2.653061224489796e-1 ...
        + xb2.*xb3.^2.*4.036999974620292e-4 ...
        - xb2.*xb3.^4.*1.946161649467695e-4 ...
        + xb2.^3.*3.993936643812048e-4 ...
        - xb2.^5.*1.946161649467695e-4
];

end