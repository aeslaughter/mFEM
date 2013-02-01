%EXAMPLE14a Beam element
% Fish and Belytschko (2009), Ex. 10.1, p. 262
%
% This is a specialized example problem, refer that replicates the example
% problem from the text exactly. For a more rubust approach refer to
% example 14b
%
% See Also EXAMPLE14b

function example14b

% Load the mFEM library
import mFEM.*

% Create the 2-element mesh
mesh = FEmesh();
mesh.addElement('Beam',[0;8]);
mesh.addElement('Beam',[8;12]);
mesh.init();

% Add flags for the boundary elements
mesh.addBoundary(1,'left');    % fixed connection
mesh.addBoundary(2,'right');   % prescribed force and rotation

% Add flag for the distributed load
mesh.addSubdomain(10, 'x<8');       % distributed load, b (integrate over element 1, using x<=8 would integrate over both elements b/c element 2 contains the node)
mesh.addSubdomain(11, 'x<=8');      % allows for dof extraction of both nodes of element 1
mesh.addSubdomain(12, 'x==8');      % allows for dof extraction of both node 1 of element 2

% Define a system
sys = System(mesh);
sys.addConstant('q',-1,'EI',10^4);
sys.addMatrix('K','EI*Ke'); 
sys.addVector('f','N''*q','Subdomain',10);

elem = mesh.element(1);
sys.addConstantVector('f',-10*elem.shape(0)','Subdomain',11, '-add');
sys.addConstantVector('f', 5, 'Subdomain', 12, 'Component', 1, '-add');
sys.addConstantVector('f', [-20,20], 'Boundary', 2, '-add');

% Create the stiffness matrix
solver = mFEM.solvers.LinearSolver(sys);
solver.addEssential('Boundary' , 1, 'Value', 0);
u = solver.solve()

% Plot the results
h = subplot(2,1,1);
%mesh.plot(u,'Axes', h, 'Component', 1, 'colorbar', 'Displacement (m)');
mesh.plot(u, 'Axes', h, 'Component', 1, '-deform', 'patch', {'EdgeColor','k'});
xlabel('x (m)');
ylabel('Displacement (m)');

h = subplot(2,1,2);
% mesh.plot(u,'Axes', h, 'Component',2, 'colorbar', 'Rotation');
mesh.plot(u, 'Axes', h, 'Component', 2, '-deform', 'patch', {'EdgeColor','k'});
xlabel('x (m)');
ylabel('Rotation (rad.)');