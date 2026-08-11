%% -----------------------------------------------------------------------
%
%  Short description: This tutorial describes how to setup a polynomial
%                     fitting using CaΣoS,
%
%
%  Problem: Given noisy data, we want to obtain a polynomial that minimizes
%           the squared error against the noisy samples. In this problem,
%           no pre-knowledge is used about the shape of the real function.
%
%           In particular, we consider samples generated from
%               
%               y = log(10*x^2+1)+w
%           
%           where w is some gaussian noise, on the domain [0,4].
%
%  License: see License file of repository
%
% -----------------------------------------------------------------------


%% Step 1) Generate/obtain noisy samples

% nominal function without noise
real_fcn = @(x) log(x.^6+1);

% obtain uniform random data from [0,4]
N_pts = 100;
pts   = 4*rand(N_pts,1);
y_w   = real_fcn(pts)+0.1*randn(N_pts,1);

%% Step 2) Setup optimization problem

% indeterminate variable
x = casos.Indeterminates('x');

% polynomial decision variable
p = casos.PS.sym('p',monomials(x,0:8));

theta = casos.PS(length(pts),1);
for i=1:length(pts)
    theta(i) = subs(p,x,pts(i))-y_w(i); 
end

% define sos problem
sos.x = p;
sos.f = theta'*theta;

% states + constraint are SOS cones
opts.Kx = struct('lin', 1);
opts.Kc = struct('sos', 0);

% build the problem
S = casos.sossol('S','mosek',sos,opts);

%% Step 3) Call solver 

% evaluate
sol = S();

% check the solution status with the unified solution status
if strcmp(S.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS')
    disp('Succesful!')
else
   disp('Unsuccesful!')
end

% extract polynomial solution
p_sol = sol.x(1);

%% Step 4) Plot the solution against the samples and the real function

% discretize the domain for evaluation
t = 0:0.1:4;
% get the polynomial function for evaluation
p_plot = p_sol.to_function;

figure(1);
hold on
% Plot samples
scatter(pts, y_w,'r.', 'DisplayName', 'Samples');
% Plot the true function without noise
plot(t, real_fcn(t), 'DisplayName', 'True function');
% Plot the polynomial approximation obtained
plot(t, full(p_plot(t)), 'DisplayName', 'Approximation without shape constraints');
legend
grid on