% A transient heat transfer example, using System class
%
% Syntax:
%   example6b
%   example6b('PropertyName', PropertyValue)
%
% Description:
%   example6b solves a simple transient heat conduction problem, with the
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

function varargout = example6b(varargin)

% Import the mFEM library
import mFEM.* mFEM.solvers.*;

% Set the default options and apply the user defined options
opt.debug = false;
opt.n = 16;
opt.element = 'Quad4';
opt = gatherUserOptions(opt,varargin{:});

% Create a FEmesh object, add the single element, and initialize it
mesh = mFEM.Mesh();
mesh.grid(opt.element,0,1,0,1,opt.n,opt.n);

% Label the boundaries
mesh.addBoundary(1); % essential boundaries (all)
% mesh.update();

% Build the system
sys = System(mesh);
sys.addConstant('D', 1 / (2*pi^2));   % thermal conductivity
sys.addMatrix('M', 'N''*N');
sys.addMatrix('K', 'B''*D*B');

% Collect the node positions for applying the essential boundary conditions
nodes = mesh.getNodes();
coord = nodes.getCoord();
x = coord(:,1);
y = coord(:,2);

% Define exact temperature
T_exact = @(x,y,t) exp(-t)*sin(pi*x).*sin(pi*y);
T = T_exact(x,y,0);
if ~opt.debug;
    mesh.plot(T); hold on;
    title('t = 0');    
end

% Plot the initial condition
if ~opt.debug;
    mesh.plot(T,'colorbar','Temperature');
    title('t = 0');
    xlabel('x');
    ylabel('y');
end

% Create and initilize the transient solver
dt = 0.1;
solver = TransientLinearSolver(sys,'dt',dt,'force',0);
solver.addEssential('Tag',1,'value',0);
solver.init(T);

% Perform 10 time-steps
for t = dt:dt:1;

    % Solve for for the temperature
    T = solver.solve();

    % Plot the results
    if ~opt.debug;
        pause(0.25);
        mesh.plot(T);
        title(['t = ', num2str(t)]);
    end
end

if opt.debug;
    varargout = {x,y,t,T};
end
