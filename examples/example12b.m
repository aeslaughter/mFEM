%EXAMPLE12b A truss example.

function example12b

% Import the mFEM library
import mFEM.*

% Create the truss structure
mesh = FEmesh();
mesh.addElement('Truss',[-1,1; 0,0]);
mesh.addElement('Truss',[0,1; 0,0]);
mesh.addElement('Truss',[1,1; 0,0]);
mesh.init();

% Label the various boundaries
mesh.addBoundary(1,{'x<=0','y>0'}); % pin connections
mesh.addBoundary(2,{'x>0','y>0'});  % roller connection
mesh.addBoundary(3,'bottom');       % applied load

% Define the known parameters
sys = System(mesh);
sys.addConstant('P',10^3,'E',10^7,'a',10^-2);
sys.addFunc('A',@(elem,x,t) elementArea(sys,elem));

% Define stiffness matrix with direct method
sys.addMatrix('K', 'A*E/L*Ke');

% Build the force vector
%vec = mFEM.kernels.ConstantVector('f', mesh, sys.get('P'), 'Boundary',3,'Component','x');
sys.addVector('f', sys.get('P'), 'Boundary',3,'Component','x');

% Create solver, apply essential boundary conditions, and solve
solver = solvers.LinearSolver(sys);
solver.addEssential({'boundary',1,'value', 0},{'boundary',2,'Component','y','value',0});
u = solver.solve();

% Plot results
mesh.plot(u,'-deform','Component','x','colorbar','x-direction disp. (m)');

function A = elementArea(sys, elem)
%ELEMENT_AREA returns the element area based on the element id

a = 10^-2;% sys.get('a');
if elem.id == 2; 
    A = 2*a; 
else
    A = a; 
end



