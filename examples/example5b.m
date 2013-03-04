% Example 8.1 from Bhatti, 2005 (p. 553)
%
% Syntax:
%   example5b
%
% Description:
%   example5b solves a simple transient heat conduction problem.
function varargout = example5b(varargin)

% Gather options
opt.debug = false;
opt = gatherUserOptions(opt, varargin{:});

% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.addElement('Tri3',1/100*[0,0; 2,0; 2,4]);
mesh.addElement('Tri3',1/100*[0,0; 2,4; 0,2]);
mesh.init();

% Label the boundaries
mesh.addBoundary(1, 'left','right');   % insulated (q = 0)
mesh.addBoundary(2, 'bottom');         % convective (q = h(T - Tinf))
mesh.addBoundary(3);                   % essential boundaries (all others)

% Create system
sys = System(mesh);

% Problem specifics
sys.addConstant('D', 3 * eye(2));      % thermal conductivity (W/(mC))
sys.addConstant('rho', 1600);          % density (kg/m^3)
sys.addConstant('c_p', 800);           % specific heat (J/kg)
sys.addConstant('h', 200);             % convective heat transfer coefficient (W/m^2)
sys.addConstant('T_s', 300);           % prescribed temperature on top (C)
sys.addConstant('T_inf', 50);          % ambient temperature (C)
sys.addConstant('T_0', 50);            % initial temperature (C)

% Add matrices
sys.addMatrix('M', 'rho*c_p*N''*N');
sys.addMatrix('K', 'B''*D*B');
sys.addMatrix('K', 'h*N''*N', 'Boundary', 2);
sys.addVector('f', 'h*T_inf*N''', 'Boundary', 2);

% Create solver
solver = solvers.TransientLinearSolver(sys, 'dt', 30);

% Add essential boundary
solver.addEssential('boundary', 3, 'value', 'T_s');

% Initialize the temperatures
T(:,1) = solver.init('T_0');

% Perform 10 time-steps
for t = 1:10;    
    % Set values for the essential boundary conditions the next time step 
    T(:, t+1) = solver.solve();
end

% Display the temperatures (in same order as p.556 of Bhatti, 2005)
if opt.debug; % debug, return M,K,f,T
    M = sys.get('M');
    K = sys.get('K');
    f = sys.get('f');
    varargout = {M.init(),K.init(),f.init(),T};
else
    ix = [1,4,2,3];
    T(ix,:)'
end
