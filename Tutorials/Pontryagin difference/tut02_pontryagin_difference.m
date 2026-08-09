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
%           a(x) = 4-x1^2-x2^2,
%           b(x) = 0.1-25*x1^2*x2^2-0.05*(x1+x2)^2.         
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
x = casos.Indeterminates('x', 2, 1);
z = casos.Indeterminates('z', 2, 1);

% define the semialgebraic sets
% a(x) >= 0 defines the set A
a = 4-x(1)^2-x(2)^2;
% b(x) >= 0 defines the set B
b = 0.1-25*x(1)^2*x(2)^2-0.05*(x(1)+x(2))^2;

% find bounding box for the set
box = find_bounding_box(a);

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
figure(1)
casos.toolboxes.sosopt.pcontour(a, 0, box, 'g')
hold on
casos.toolboxes.sosopt.pcontour(b, 0, box, 'b')
casos.toolboxes.sosopt.pcontour(c_sol, 0, box, 'r')
legend('a(x)=0', 'b(x)=0', 'c(x)=0')
grid on