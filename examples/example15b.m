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
opt.debug = false;
opt = gatherUserOptions(opt,varargin{:});

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Element',opt.element,'-time');
mesh.grid(0,1,0,1,0,1,opt.n,opt.n,opt.n);
mesh.init();

% Label the boundaries
mesh.addBoundary(1); % essential boundaries (all)

% Build the system
sys = System(mesh);
sys.addConstant('D', 1 / (2*pi^2), 'dt', 0.1); 
sys.addMatrix('M', 'N''*N');
sys.addMatrix('K', 'B''*D*B');

% Define exact temperature
[x,y,z] = mesh.getNodes();
T_exact = @(x,y,z,t) exp(-t)*sin(pi*x).*sin(pi*y).*sin(pi*z);
T = T_exact(x,y,z,0);

% Plot the initial condition
if ~opt.debug;
    options = {'slice',{'x',[0.5,1],'y',[0.5,1],'z',[0,0.2]},...
               'colorbar','Temperature','-ShowElements'};
    mesh.plot(T,'-new',options{:});
    title('t = 0'); xlabel('x'); ylabel('y');
end

% Create and initilize the transient solver
solver = TransientLinearSolver(sys, 'dt', sys.get('dt'), 'force', 0,'-disableAll');
solver.addEssential('Boundary', 1, 'value', 0);
solver.init(T);

% Perform 10 time-steps
for t = 1:10;
    T = solver.solve();
    if ~opt.debug;
        pause(0.25);
        mesh.plot(T, options{:});
        title(['t = ', num2str(sys.get('dt')*t)]);
    end
end

max(abs(T - T_exact(x,y,z,sys.get('dt')*t)))