function box = find_bounding_box(a, varargin)
% Description: Find axis-aligned bounding box of polynomial a(x) returning 
%              the bounding box [x1_min, x1_max, x2_min, x2_max, ...] for 
%              the polynomial a(x).

% process input options
if isempty(varargin)
    param = struct();
else
    param = varargin{:};
end

if ~isfield(param,'deg_s'), param.deg_s = 1; end
if ~isfield(param,'sdp_solver'), param.sdp_solver = 'mosek'; end

% get indeterminates and dimension
x = a.indeterminates;
n = length(x); 

% pre-allocate bounding box
box = zeros(1, 2*n);

% define auxiliary variables
lvl = casos.PS.sym('lvl', 1);
sgn = casos.PS.sym('sgn', 1);
cxv = casos.PS.sym('cxv', n);

% sos multiplier
s = casos.PS.sym('s', monomials(x,0:param.deg_s), 'gram');

% build sos problem
sos.f = sgn*lvl;
sos.g = -s*a+sgn*(lvl-cxv'*casos.PS(x));
sos.x = [lvl; s];
sos.p = [sgn; cxv];

% set cones
opts = struct('Kx', struct('lin', 1, 'sos', 1), ...
              'Kc', struct('sos', 1));

% build solver
S = casos.sossol('S', param.sdp_solver, sos, opts);          

% solve for all dimensions and signs
for i=1:n
    tec = zeros(n, 1);
    tec(i) = 1;
    sd = 2*(i-1);
    for k=[-1, 1] 
        sol = S('p', [k;tec]);

        % check status
        if ~strcmp( S.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
            error('Solver failed while determing the bounding box!')
        end

        % save bound
        id = 0.5*k+1.5;
        box(sd+id) = full(sol.x(1));
    end
end

end