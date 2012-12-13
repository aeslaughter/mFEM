% A transient heat transfer example, using System class
%
% Syntax:
%   example6c
%   example6c('PropertyName', PropertyValue)
%
% Description:
%   example6c solves a simple transient heat conduction problem, with the
%   default settings.
%
%   example6c('PropertyName', PropertyValue) allows the user to customize
%   the behavior of this example using the property pairs listed below.
%
% Example6c Property Descriptions
%
% N
%   scalar
%    The number of elements in the x and y directions, the default is 32.
%
% Element
%   {'Quad4'} | 'Tri3' | 'Tri6'
%   Specifies the type of element for the mesh

function example6c(varargin)

% Import the mFEM library
import mFEM.* mFEM.solvers.*;

% Set the default options and apply the user defined options
opt.n = 16;
opt.element = 'Quad4';
opt = gather_user_options(opt,varargin{:});

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Element',opt.element);
mesh.grid(0,1,0,1,opt.n,opt.n);
mesh.init();

% Label the boundaries
mesh.add_boundary(1); % essential boundaries (all)

% Build the system
sys = System(mesh);
sys.add_constant('D', 1 / (2*pi^2));   % thermal conductivity
sys.add_matrix('M', 'N''*N');
sys.add_matrix('K', 'B''*D*B');

% Collect the node positions for applying the essential boundary conditions
nodes = mesh.get_nodes();
x = nodes(:,1);
y = nodes(:,2);

% Define exact temperature
T_exact = @(x,y,t) exp(-t)*sin(pi*x).*sin(pi*y);
T = T_exact(x,y,0);

% Plot the initial condition
figure; hold on;
mesh.plot(T);
title('t = 0');
xlabel('x');
ylabel('y');
cbar = colorbar;
set(get(cbar,'YLabel'),'String','Temperature');

% Create and initilize the transient solver
dt = 0.1;
solver = TransientLinearSolver(sys, 'dt', dt, 'force', 0);
solver.add_essential_boundary('id', 1, 'value',0);
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
