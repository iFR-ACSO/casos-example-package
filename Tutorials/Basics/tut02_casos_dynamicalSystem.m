%% -----------------------------------------------------------------------
%
% Short description: This tutorial describes how to use CaΣoS to setup a
%                    simple nonlinear dynamical system, compute a linear 
%                    control law and setup a (simple) simulation.
%
%   License: see License file of repository
%
% -----------------------------------------------------------------------


%% Step 1) Define open-loop nonlinear system
% indeterminate variable
x = casos.Indeterminates('x');
u = casos.Indeterminates('u');

% polynomial dynamics
f = x^4 + 10*x + u;

%% Step 2) Compute an LQR controller

% Step 2.1: linearize system using CaSoS
A = nabla(f,x);
B = nabla(f,u);

A0 = full(subs(A,[x;u],zeros(2,1)));
B0 = full(subs(B,[x;u],zeros(2,1)));

% Step 2.2 let us design a simple LQR controller
[K,P] = lqr(A0,B0,1,1);


% Step 2.3 generate the closed-loop dynamics
f_cl = subs(f,u,-K*x);

% Step 2.4 make a casos function of the closed-loop system for repeated evaluation
f_cl_fun = to_function(f_cl);

%% Step 3) simulation (simple)
dt    = 0.1; % step size
T     = 10; % simulation duration
steps = linspace(0,T,T/dt);
% allocate
x     = zeros(1,length(steps));

% inital state
x(1) = 1;

for k = 2:length(steps)
 % Euler-forward
 x(k) = x(k-1) + dt*full(f_cl_fun(x(k-1)));

end

% plotting
figure('Name','Closed-loop Trajectory')
plot(steps,x)


