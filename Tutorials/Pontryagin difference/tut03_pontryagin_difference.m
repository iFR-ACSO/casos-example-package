% -------------------------------------------------------------------------
%
%  Short description: This tutorial describes how to approximate the 
%                     Pontryagin difference between semialgebraic sets
%                     using a sum-of-squares-based procedure.
%
%  Problem: Find a polynomial c(x) that approximates the Pontryagin 
%           difference between sets defined by a(x) >= 0 and b(x) >= 0.
%           The optimization problem is:
%
%           minimize ∫ c(x) dx
%           subject to:
%           a(x+z) - c(x) - s(x,z)*b(z) >= 0  (SOS constraint)
%           where s(x,z) is an SOS multiplier
%
%           c(x) is a polynomial decision variable 
%           s(x,z) is an SOS multiplier 
%   
%           The example is based on and taken from [1], where
%           
%           a(x) = 1-x1^6-x2^6-x3^6+5*x1^4*x2*x3-3*x1^4*x2^2+
%                  -10*x1^2*x2^3*x3-3*x1^2*x2^4+x2^5*x3,
%           b(x) = 0.0001-x1^6-x2^6-x3^6.         
%
%  Reference
%     [1] A. Cotorruelo, I. Kolmanovsky, and E. Garone, "A sum-of-squares-based 
%         procedure to approximate the Pontryagin difference of basic 
%         semi-algebraic sets", in Automatica, vol. 135, pp. 109783, 2022.
%
%  License: see License file of repository
%
% -------------------------------------------------------------------------

%% Define indeterminates and problem parameters
% indeterminate variable
x = casos.Indeterminates('x', 3, 1);
z = casos.Indeterminates('z', 3, 1);

% define the semialgebraic sets
% a(x) >= 0 defines the set A
a = 1-x(1)^6-x(2)^6-x(3)^6+5*x(1)^4*x(2)*x(3)-3*x(1)^4*x(2)^2-10*x(1)^2*x(2)^3*x(3)-3*x(1)^2*x(2)^4+x(2)^5*x(3);
% b(x) >= 0 defines the set B
b = 1e-4-x(1)^6-x(2)^6-x(3)^6;

% find bounding box for the set
param.deg_s = 2;
box = find_bounding_box(a, param);

%% Define decision variables
% polynomial c(x) - approximates the Pontryagin difference
c = casos.PS.sym('c', monomials(x, 0:10)); 

% SOS multiplier s(x,z)
s = casos.PS.sym('s', monomials([x;z], 0:2), 'gram'); 

%% Setup the problem struct
% Pontryagin difference constraint: a(x+z)-c(x)-s(x,z)*b(z)>= 0
P = subs(a,x,x+z)-c-s*subs(b,x,z);

% cost function: minimize ∫ c(x) dx over the bounding box
cost = int(c,x(1), box(1:2));
cost = int(cost, x(2), box(3:4));
cost = int(cost, x(3), box(5:6));

% setup the problem struct
sos = struct('x', [c;s], 'f', -cost, 'g', P);

% provide the problem size, i.e., size of cones
opts.Kx = struct('lin', 1, 'sos', 1);
opts.Kc = struct('sos', 1);

% solver options for underlying convex SOS problem
% if true, error returns infeasible
opts.error_on_fail = false; 
% turn off Newton polytope reduction
opts.newton_solver = 'mosek';    

% choose the underlying SDP solver
solver = 'mosek';

% generate a CaΣoS solver instance
S = casos.sossol('S', solver, sos, opts);

%% Call the solver and extract solution
% call solver
startTimeSolver = tic;
sol = S();
elapsedTimeSolver = toc(startTimeSolver);
fprintf('Solver: elapsed time %f seconds\n', elapsedTimeSolver)

% check the solution status with the unified solution status
if strcmp(S.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
    fprintf('Optimal solution found\n')
else
    disp('Unsuccesful!')
end

% extract the solution
c_sol = sol.x(1);
s_sol = sol.x(2);

%% Visualize the results
% plot surfaces
fig = figure(1);
casos.toolboxes.sosopt.pcontour3(a, 0, box);
hold on
casos.toolboxes.sosopt.pcontour3(b, 0, box);
casos.toolboxes.sosopt.pcontour3(c_sol, 0, box);

% change settings
grh = fig.Children.Children;
sg  = findobj(grh, 'Type', 'patch');
sg(1).FaceAlpha = 0.5;  sg(1).FaceColor = 'red';
sg(2).FaceAlpha = 0.5;  sg(2).FaceColor = 'blue';
sg(3).FaceAlpha = 0.25; sg(3).FaceColor = 'green';
legend('a(x)=0', 'b(x)=0', 'c(x)=0')
grid on