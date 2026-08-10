%% -----------------------------------------------------------------------
%
% Example 2: Stabilization of a cubesat in circular orbit
%
% Short description: Verify almost global stability of a cubesat in a
%                    circular orbit subject to gravity gradient torques
%                    using quaternions for attitude parametrization.
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
xbar       = casos.PD('xb', 7);
xbar_star1 = [1; zeros(6,1)];     % stable equilibrium in original coordinates
xbar_star2 = [-1; zeros(6,1)];    % unstable equilibrium in original coordinates

% xbar(1:4): attitude quaternion of B w.r.t. O (scalar part first)
% xbar(5:7): rotational rates expressed in body coordinates

% Gravity gradient torque and the quaternion feedback regulator are already
% contained in the generated dynamics.
fbar = example2_dynamics(xbar);

% Manifold constraint in original coordinates
hbar = xbar(1:4)'*xbar(1:4) - 1;


%% Step 2) Shift and scale system dynamics and define SOS quantities

% Coordinate transformation:
%
%   x = S*(xbar - xbar_star1)
%
S = diag([ones(4,1); 15.*ones(3,1)]);

% Transformed state
x = casos.PD('x', 7);

% Inverse map:
%
%   xbar = S\x + xbar_star1
%
% used to substitute into the original expressions.
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


% Enforce positive definiteness on V
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
    1e-3*max(sol.x(2).poly2basis));


%% Step 8) Create plot

if ~plot_flag
    return
end

% Lyapunov function in original coordinates
Vbar_sol = subs(Vsol, x, S*(xbar - xbar_star1));


% Random initial condition
%
% x0 = -1 + 2.*rand(7,1);
% x0(1:4) = [-.93; 0.1; 0; 0];
% x0(1:4) = x0(1:4)/norm(x0(1:4)); % quaternion unit constraint
% x0(5:7) = 0.05.*x0(5:7);         % scale rates to reasonable value


% Initial condition from paper plot
x0 = [
    -0.9943
     0.1069
     0
     0
     0.0377
     0.0043
     0.0363
];


% Simulate the true, untransformed dynamics
[t,y] = ode45( ...
    @(~,xx) example2_dynamics(xx), ...
    linspace(0,250,701), ...
    x0);

t = t.';
y = y.'; % states in rows


% Evaluate Lyapunov function
Vbar_solFun = Vbar_sol.to_function();

yrows = num2cell(y,2);

Vbar_val = full(Vbar_solFun(yrows{:}));


% Plot states and normalized Lyapunov function
figure('Name', 'States and Lyapunov function (original coordinates)');
clf;
hold on;

plot( ...
    t, ...
    y(1:4,:), ...
    LineStyle='--', ...
    LineWidth=.8);

plot( ...
    t, ...
    rad2deg(y(5:7,:)), ...
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
    '$\bar{x}_4$', ...
    '$\bar{x}_5 \,[^{\circ}/s]$', ...
    '$\bar{x}_6 \,[^{\circ}/s]$', ...
    '$\bar{x}_7 \,[^{\circ}/s]$', ...
    '$\bar{V}_{\mathrm{s}}(\bar{x})$', ...
    'Interpreter', 'latex', ...
    'Location', 'eastoutside');

set( ...
    gca, ...
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 12);

hold off;


function f = example2_dynamics(in1)
%EXAMPLE2_DYNAMICS
%    F = EXAMPLE2_DYNAMICS(IN1)
%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    04-Aug-2026 09:28:14
%Gravity-gradient stabilisation of a cubesat in circular orbit (Example 2).
%
%States:
%x(1:4)  attitude quaternion of B w.r.t. O (scalar part first)
%x(5:7)  body rates [rad/s]
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
xb7 = in1(7,:);

t2 = xb1.^2;
t3 = xb1.^3;
t4 = xb2.^2;
t5 = xb3.^2;
t6 = xb4.^2;

mt1 = [
    xb2.*xb5.*-5.0e-1 ...
        - xb3.*xb6.*5.0e-1 ...
        - xb4.*xb7.*5.0e-1;

    xb1.*xb5.*5.0e-1 ...
        + xb3.*xb7.*5.0e-1 ...
        - xb4.*xb6.*5.0e-1;

    xb1.*xb6.*5.0e-1 ...
        - xb2.*xb7.*5.0e-1 ...
        + xb4.*xb5.*5.0e-1;

    xb1.*xb7.*5.0e-1 ...
        + xb2.*xb6.*5.0e-1 ...
        - xb3.*xb5.*5.0e-1;

    xb2.*-6.4e-3 ...
        - xb5.*8.0e-2 ...
        - xb7.*1.106816514833168e-3 ...
        + t2.*xb7.*2.213633029666336e-3 ...
        + t5.*xb7.*2.213633029666336e-3 ...
        + xb1.*xb2.*xb6.*2.213633029666336e-3 ...
        - xb3.*xb4.*xb6.*2.213633029666336e-3;

    xb3.*-6.4e-3 ...
        - xb6.*8.0e-2 ...
        + t3.*xb3.*3.523742373754575e-3 ...
        - xb1.*xb3.*1.761871186877288e-3 ...
        + xb2.*xb4.*1.761871186877288e-3 ...
        + xb5.*xb7.*2.653061224489796e-1 ...
        - xb2.*xb4.^3.*3.523742373754575e-3 ...
        - t2.*xb2.*xb4.*3.525042419172338e-3 ...
        - t4.*xb1.*xb3.*1.300045417762998e-6 ...
        + t6.*xb1.*xb3.*3.525042419172338e-3 ...
        + t5.*xb2.*xb4.*1.300045417762998e-6 ...
        - xb1.*xb2.*xb5.*1.626342634040573e-3 ...
        - xb1.*xb4.*xb7.*2.800923425292098e-3 ...
        - xb2.*xb3.*xb7.*2.800923425292098e-3 ...
        + xb3.*xb4.*xb5.*1.626342634040573e-3
];

mt2 = [
    xb4.*-6.4e-3 ...
        + xb5.*8.131713170202867e-4 ...
        - xb7.*8.0e-2 ...
        - t2.*xb5.*1.626342634040573e-3 ...
        - t3.*xb4.*1.300045417762998e-6 ...
        - t5.*xb5.*1.626342634040573e-3 ...
        + xb1.*xb4.*6.500227088814988e-7 ...
        + xb2.*xb3.*6.500227088814988e-7 ...
        - xb5.*xb6.*2.653061224489796e-1 ...
        - xb2.*xb3.^3.*1.300045417762998e-6 ...
        - t2.*xb2.*xb3.*3.525042419172338e-3 ...
        + t4.*xb1.*xb4.*3.523742373754575e-3 ...
        - t5.*xb1.*xb4.*3.525042419172338e-3 ...
        + t6.*xb2.*xb3.*3.523742373754575e-3 ...
        + xb1.*xb4.*xb6.*2.800923425292098e-3 ...
        + xb2.*xb3.*xb6.*2.800923425292098e-3
];

f = [mt1; mt2];

end