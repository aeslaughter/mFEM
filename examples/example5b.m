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

% Create a Mesh object, add the single element, and initialize it
mesh = mFEM.Mesh();
mesh.createNode(1/100*[0,0; 0,2; 2,0; 2,4]);
mesh.createElement('Tri3', [1,3,4; 1,4,2]);
mesh.init();

% Label the boundaries
mesh.addBoundary('insulated','left','right');  
mesh.addBoundary('convective','bottom'); 
mesh.addBoundary('essential','y>=0.02');
mesh.update();

% Create system
sys = mFEM.System(mesh);

% Problem specifics
sys.addConstant('D', 3 * eye(2));      % thermal conductivity (W/(mC))
sys.addConstant('rho', 1600);          % density (kg/m^3)
sys.addConstant('c_p', 800);           % specific heat (J/kg)
sys.addConstant('h', 200);             % convective heat transfer coefficient (W/m^2)
sys.addConstant('T_s', 300);           % prescribed temperature on top (C)
sys.addConstant('T_inf', 50);          % ambient temperature (C)
sys.addConstant('T_0', 50);            % initial temperature (C)

% Add matrices
sys.addMatrix('M','rho*c_p*N''*N');
sys.addMatrix('K','B''*D*B');
sys.addMatrix('K','h*N''*N','Tag','convective');
sys.addVector('f','h*T_inf*N''','Tag','convective');

% Create solver
solver = mFEM.solvers.TransientLinearSolver(sys,'dt',30);

% Add essential boundary
solver.addEssential('Tag','essential','Value','T_s');

% Initialize the temperatures
T(:,1) = solver.init('T_0');

% Perform 10 time-steps
for t = 1:10;    
    % Set values for the essential boundary conditions the next time step 
    T(:, t+1) = solver.solve();
end

% Display the temperatures (in same order as p.556 of Bhatti, 2005)
if opt.debug; % debug, return M,K,f,T
    M = sys.assemble('M'); 
    K = sys.assemble('K'); 
    f = sys.assemble('f'); 
    varargout = {M,K,f,T};
else
    T'
end
