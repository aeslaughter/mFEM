% A transient heat transfer example, using System class
%
% Syntax:
%   example15b
%   example15b('PropertyName', PropertyValue)
%
% Description:
%   example15b solves a simple transient heat conduction problem, with the
%   default settings.
%
%   example15c('PropertyName', PropertyValue) allows the user to customize
%   the behavior of this example using the property pairs listed below.
%
% Example15c Property Descriptions
%
% N
%   scalar
%    The number of elements in the x,y, and z directions, the default is 8.
%
% Element
%   {'Hex8'}
%   Specifies the type of element for the mesh

function example15b(varargin)

% Import the mFEM library
import mFEM.* mFEM.solvers.*;

% Set the default options and apply the user defined options
opt.n = 8;
opt.element = 'Hex8';
opt = gatherUserOptions(opt,varargin{:});

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Element',opt.element);
mesh.grid(0,1,0,1,0,1,opt.n,opt.n,opt.n);
mesh.init();

% Label the boundaries
mesh.addBoundary(1); % essential boundaries (all)

% Build the system
sys = System(mesh);
sys.addConstant('D', 1 / (2*pi^2));   % thermal conductivity
sys.addMatrix('M', 'N''*N');
sys.addMatrix('K', 'B''*D*B');

% Collect the node positions for applying the essential boundary conditions
nodes = mesh.getNodes();
x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);

% Define exact temperature
T_exact = @(x,y,z,t) exp(-t)*sin(pi*x).*sin(pi*y).*sin(pi*z);
T = T_exact(x,y,z,0);

% Plot the initial condition
mesh.plot(T,'-new');
title('t = 0');
xlabel('x');
ylabel('y');

% Create and initilize the transient solver
dt = 0.1;
solver = TransientLinearSolver(sys, 'dt', 0.1, 'force', 0);
solver.addEssential('Boundary', 1, 'value',0);
solver.init(T);

% Perform 10 time-steps
for t = dt:dt:1;

    % Solve for for the temperature
    T = solver.solve();

    % Plot the results
    pause(0.25);
    mesh.plot(T);
    title(['t = ', num2str(t)]);
end
