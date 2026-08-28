%% -----------------------------------------------------------------------
%
% Short description: This tutorial shows how to plot polynomial level sets.
%                    It demonstrates how to set it up using legacy code
%                    from sosopt and using fcontour. It demonstrates how to
%                    plot 2D and 3D example and how to plot slices for
%                    higher dimensional systems.
% 

%
%   License: see License file of repository
%
% -----------------------------------------------------------------------

% import legacy functions
import casos.toolboxes.sosopt.pcontour
import casos.toolboxes.sosopt.pcontour3


%% 2D example
% Indeterminate variable 
y = casos.Indeterminates('y',2);


% assume a simple ball
p1 = y(1)^2 + y(2)^2 - 1;

% using legacy code
figure('Name','Plot level set with sosopt legacy code of pcontour')
pcontour(p1, ...            % polynomials
         0, ...             % level set; here zero sublevel set
         [-2 2 -2 2],...    % domain (min/max) for each dimension
         'r--' )            % set color and line type
legend('Zero level set')



% using fcontour
figure('Name','Manually plot the 2D level set')
% generate a function
pfun = to_function(p1);

fcontour(@(x,y) ...
        full(pfun(x,y)), ... % evaluate polynomial (full() to get double)
        [-2 2], ...          % domain for plotting
        'b--', ...           % set color and line type
        "LevelList", [0 0])  % set level set
legend('Zero level set')
xlabel('x_1')
ylabel('x_2')


%% 3D example
% Indeterminate variable 
y = casos.Indeterminates('y',3);


% assume a simple ball
p2 = y(1)^2 + y(2)^2  + y(3)^2 - 1;

% using legacy code
figure('Name','Plot 3D level set with sosopt legacy code of pcontour3')
pcontour3(p2, ...            % polynomials
    0, ...                   % level set; here zero sublevel set
    [-2 2 -2 2 -2 2])        % domain 3D
legend('Zero level set')


% using fcontour
figure('Name','Manually plot the 3D level set')
% generate a function
pfun3 = to_function(p2);

x = casadi.SX.sym('x',p2.nvars);

xcell = num2cell(x);
W_expr_old = pfun3(xcell{:});
W_expr = casadi.substitute(W_expr_old, x, x);
W_fun = casadi.Function('W_fun',{x},{W_expr});

fimplicit3(@(x,y,z) ...
           full(pfun3(x,y,z)), ...  % evaluate polynomial (full() to get double)  
            [-2 2 -2 2 -2 2], ...   % domain 3D            
            'EdgeColor','none', ... % Define appearance
             'FaceAlpha',0.6)
colormap(parula)
camlight headlight
lighting gouraud
axis equal
legend('Zero level set')
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')

%% 5D example
% Indeterminate variable 
y = casos.Indeterminates('y',5);

% assume a polynomial in 5 indeterminates
p3 = 2*y(1)^2 + y(2)^2 + 3*y(3)^2 + y(4)^2 +  1.5*y(5)^2    - 1;

% ------------------------------------------------------------
% p3 cannot be plotted directly, we need slices or projections
% ------------------------------------------------------------

% Plot slice of x1-x2:
% substitute fixed values for x3,x4,x5; here all are set to zero
p3_x1x2_slice = subs(p3,y(3:end),zeros(3,1));

% ------------------------------------------------------------
% Note: 
%   - For a slice the values must not be zero!
%   - To plot e.g. the x4-x5 level set set the other 
%     indeterminates to a desired value (e.g. all to zero. 
% ------------------------------------------------------------

% using legacy code to plot the slice x1-x2
figure('Name','Plot slice using legacy code')
pcontour(p3_x1x2_slice , ...            % polynomials
    0, ...             % level set; here zero sublevel set
    [-2 2 -2 2],...    % domain (min/max) for each dimension
    'r--' )            % set color and line type
legend('Zero level set')



% using fcontour
figure('Name','lot slice anually plot the 2D level set')
% generate a function
pfun = to_function(p3_x1x2_slice );

fcontour(@(x,y) ...
    full(pfun(x,y)), ... % evaluate polynomial (full() to get double)
    [-2 2], ...          % domain for plotting
    'b--', ...           % set color and line type
    "LevelList", [0 0])  % set level set
legend('Zero level set')
xlabel('x_1')
ylabel('x_2')