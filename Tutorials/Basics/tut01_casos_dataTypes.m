%% -----------------------------------------------------------------------
%
% Short description: This tutorial describes the different CaΣoS data
%                    types, how to setup them up and how to mainpulate
%                    them.
%
% 
%   Indeterminate variables             casos.Indeterminates
%   Symbolic polynomial variables       casos.PS.sym
%   Nmerical polynomial variables       casos.PD
%
%   License: see License file of repository
%
% -----------------------------------------------------------------------


%% Indeterminate variable 
% Are used, for example, to setup the system dynamics
% scalar
x = casos.Indeterminates('x');
% vector
y = casos.Indeterminates('y',2);

% polynomial dynamics 
f1 = x^4 + 10*x;

f2 = y(1)^2*y(2) + y(1);

%% Symbolic variables and polynomials

% scalar, for example, level set
g = casos.PS.sym('g');

% a univariate polynomial up to degree 2
p = casos.PS.sym('p',monomials(x,0:2));

% a univariate polynomial with custom monomials
p1 = casos.PS.sym('p',monomials([x x^2 x^4]));

% multivariate polynomial
p2 = casos.PS.sym('p',monomials(y,0:2));

% multivariate polynomial with custom monomials
p3 = casos.PS.sym('p',monomials([y(1)*y(2) y(1)^2 y(2)^4]));

% n = 2 multivariate polynomials
p4 = casos.PS.sym('p',monomials(y,0:2), [2 1]);

% n x n with n = 2 multivariate polynomial matrix
p5 = casos.PS.sym('p',monomials(y,0:2), [2 2]);

%% getting the coefficients and sparsity pattern of a polynomial

[coeffs_p2,monomials_p2] = poly2basis(p2);

% from symbolic to numerical values
p2_num = casos.PD(monomials_p2,rand(monomials_p2.nnz,1));

% get the numerical coefficients only
[coeffs_p2__num,~] = poly2basis(p2_num);

%% manipulating numerical polynomials

% remove small coefficients
p2_num = casos.PD(monomials_p2,(1 - 10^-2)*randn(monomials_p2.nnz,1)+(10^-2))

% remove all coefficients smaller than 1
p2_num_reduced = remove_coeffs(p2_num,1)