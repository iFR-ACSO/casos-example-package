%--------------------------------------------------------------------------
%
% Short Description:   This script implements the synthesis of an
% input-to-state (ISS) Lyapunov function. The example is based on the method from
% [1], but extended to also synthesis a (liner) control law.
%
%  Date:    21.03.2026
%
%
%  Reference:
%  [1]	H. Ichihara, „Sum of Squares Based Input-to-State Stability Analysis
%       of Polynomial Nonlinear Systems", SICE Journal of Control, Measurement,
%       and System Integration, Bd. 5, Nr. 4, S. 218–225, July 2012,
%       doi: 10.9746/jcmsi.5.218.
%
%   License: see License file of repository
%
%------------------------------------------------------------------

clc
clear

%% problem definition
% system states
x = casos.PS('x',6);
u = casos.PS('u',3);
w = casos.PS('w',3);
s = casos.PS('s',1);

%% satellite dynamics (simple)

% simple bounds on rates;
omegaMax1 = 0.5*pi/180;
omegaMax2 = 0.2*pi/180;
omegaMax3 = 0.2*pi/180;

x_low =  [-omegaMax1 -omegaMax2 -omegaMax3]';
x_up  =  [ omegaMax1  omegaMax2  omegaMax3]';

% scaling matrix for system states
Dx   = diag([1/(x_up(1)-x_low(1)),1/(x_up(2)-x_low(2)),1/(x_up(3)-x_low(3)),1,1,1]);

Dxin = inv(Dx);

% inertia tensor (simple)
J = diag([1;1;1]);

% cross-product matrix
cpm = @(x) [   0  -x(3)  x(2);
    x(3)   0   -x(1);
    -x(2)  x(1)   0 ];

% MRP derivative
B = @(sigma) (1-sigma'*sigma)*eye(3)+ 2*cpm(sigma)+ 2*sigma*sigma';

% dynamics
f =  [-J\cpm(x(1:3))*J*x(1:3) + J\u; % omega_dot
    1/4*B(x(4:6))*x(1:3)];         % sigma_dot


% ISS Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2:4));
K = casos.PS.sym('k',monomials(x,1),[3 1]); 

% closed-loop dynamics
f = subs(f,u,K);

% scale dynamics
f = Dx*subs(f,x,inv(Dx)*x);


% K_inf (see Example 1)
a      = casos.PS.sym('ca',monomials([s^2 s^4]));
a_lbar = casos.PS.sym('cau',monomials([s^2 s^4]));
a_bar  = casos.PS.sym('cao',monomials([s^2 s^4]));
sigma  = casos.PS.sym('si',monomials([s^2 s^4]));


% constraints
[c_a_lbar,~]  = poly2basis(a_lbar); % helper to get "norm"
c_a_lbar      = casos.PS(c_a_lbar);

% lower bound on ISS-Lyapunov
g1 = V - c_a_lbar(1)*(x'*x) - c_a_lbar(2)*(x'*x)^2;

[c_a_bar,~]  = poly2basis(a_bar); % helper to get "norm"
c_a_bar      = casos.PS(c_a_bar);


% upper bound on ISS-Lyapunov
g2 = c_a_bar(1)*(x'*x) + c_a_bar(2)*(x'*x)^2 - V;

[c_sigma ,~] = poly2basis(sigma ); % helper to get "norm"
c_sigma      = casos.PS(c_sigma);  

[c_a,~]      = poly2basis(a);      % helper to get "norm"
c_a          = casos.PS(c_a);


% ISS-Lyapunov dissipation inequality: −∇𝑉(𝑥)𝑓(𝑥,𝑤) + 𝜎(‖𝑤‖) − 𝛼(‖𝑥‖) is SOS
g3 = c_sigma(1)*(w'*w) + c_sigma(2)*(w'*w)^2 - nabla(V,x)*f - c_a(1)*(x'*x) - c_a(2)*(x'*x)^2  ;

% ensure univariate polynomials are Κ∞ functions 
g4 = s*nabla(a_bar,s);
g5 = s*nabla(a_lbar,s);
g6 = s*nabla(sigma,s);
g7 = s*nabla(a,s);

% cost function
f = dot(V,V);

% parameter
p = [];

% polynomial decision variables
x_lin   = [V; K;a_bar; a_lbar; a; sigma]; % linear decision variables
x_sos   = [];                             % sos decision variables


% setup constraints
g_lin = [];                      % linear constraints
g_sos =  [g1;g2;g3;g4;g5;g6;g7]; % SOS constraints


% setup SOS problem struct
if ~isempty(f)
    sos = struct('x',[x_lin;x_sos], ...
        'f', f,...
        'p',p,...
        'g',[g_lin;g_sos]);

else
    sos = struct('x',[x_lin;x_sos], ...
        'p',p,...
        'g',[g_lin;g_sos]);
end

% Provide the problem size i.e. size of cones
nx_sos = length(x_sos);
nx_lin = length(x_lin);

ng_sos = length(g_sos);
ng_lin = length(g_lin);

opts.Kx.sos = nx_sos;
opts.Kx.lin = nx_lin;

opts.Kc.lin = ng_lin;
opts.Kc.sos = ng_sos;

% options
opts.sossol_options.newton_solver = [];


%% get solver
sdp_solver = 'mosek'; % 'mosek', 'scs', 'clarabel'


S = casos.nlsossol('S1','sequential',sos,opts);

%% solve
x0 = casos.PD([x'*x; ones(3,1)*(x'*x);1;1;1;1]);

sol = S('x0',x0);

Vsol = remove_coeffs(subs(sol.x(1),x,Dx*x),1e-6)
Ksol = remove_coeffs(subs(sol.x(2:4),x,Dx*x),1e-6)
Ksol = remove_coeffs(subs(sol.x(2:4),x,Dx*x),1e-6);
