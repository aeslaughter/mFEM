%EXAMPLE14c Beam element
% Fish and Belytschko (2009), Ex. 10.1, p. 262
%
% See Also EXAMPLE14c

function example14c(varargin)

% User options
opt.n = 10;
opt.direct = false;
opt.debug = false;
opt = gatherUserOptions(opt, varargin{:});

% Load the mFEM library
import mFEM.*

% Create the 2-element mesh
mesh = FEmesh('Element','Beam');
mesh.grid(0,12,opt.n*3); % using ensures that nodes always are at nodes
mesh.init();

% Add flags for the boundary elements
mesh.addBoundary(1, 'left');    % fixed connection
mesh.addBoundary(2, 'right');   % prescribed force and moment

% Add flag for the distributed load and point loads
mesh.addSubdomain(10, 'x<8');    % distributed load, b
mesh.addSubdomain(11, 'x==4');   % point load, P1
mesh.addSubdomain(12, 'x==8');   % point load, P2

% Define a system
sys = System(mesh);
sys.addConstant('b',-1,'EI',10^4);
sys.addVector('f','N''*b','Subdomain',10);

if opt.direct;
    sys.addMatrix('K','EI*Ke','-direct'); 
else
    sys.addMatrix('K','EI*B''*B');
end

sys.addConstantVector('f', -10,'Subdomain', 11, 'Component', 1);
sys.addConstantVector('f', 5, 'Subdomain', 12, 'Component', 1, '-add');
sys.addConstantVector('f', [-20,20], 'Boundary', 2, '-add');

% Create the stiffness matrix
solver = mFEM.solvers.LinearSolver(sys);
solver.addEssential('Boundary' , 1, 'Value', 0);
u = solver.solve();

% Plot the results
h = subplot(4,1,1);
mesh.plot(u,'Axes', h, 'Component', 1, '-deform', 'patch', {'EdgeColor','k','LineWidth',2});
xlabel('x (m)');
ylabel('Displacement (m)');

h = subplot(4,1,2);
mesh.plot(u, 'Axes', h, 'Component', 2, '-deform', 'patch', {'EdgeColor','k','LineWidth',2});
xlabel('x (m)');
ylabel('Rotation (rad.)');

% Compute the moment and shear
s = []; sx = []; m = []; mx = [];
for e = 1:mesh.n_elements;
    
    % Current element
    elem = mesh.element(e);
    dof = elem.getDof();  
    
    % Create functions for N and B
    B = @(xi) elem.shapeDeriv(xi);
    D = elem.dNdx3();

    % Shear
    s(end+1) = sys.get('EI')*D*u(dof);
    s(end+1) = s(end);
    sx(end+1) = elem.nodes(1);
    sx(end+1) = elem.nodes(2);
    
    % Moments at quass points
    for i = 1:length(elem.qp);
        m(end+1) = sys.get('EI')*B(elem.qp{i})*u(dof);
        mx(end+1) = elem.getPosition(elem.qp{i});
    end
end

% Plot the shear force
subplot(4,1,3);
plot(sx,s,'b-', 'LineWidth', 2);
xlabel('x (m)');
ylabel('Shear (N)');

% Plot the moment
subplot(4,1,4);
plot(mx,m,'b-', 'LineWidth', 2);
xlabel('x (m)');
ylabel('Moment (Nm)');
