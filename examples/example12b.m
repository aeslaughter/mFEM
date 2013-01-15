%EXAMPLE12b A truss example.

function example12b

% Import the mFEM library
import mFEM.*

% Create the truss structure
mesh = FEmesh('Element','Truss');
mesh.add_element([-1,1; 0,0]);
mesh.add_element([0,1; 0,0]);
mesh.add_element([1,1; 0,0]);
mesh.init();

% Label the various boundaries
mesh.add_boundary(1,{'x<=0','y>0'}); % pin connections
mesh.add_boundary(2,{'x>0','y>0'});  % roller connection
mesh.add_boundary(3,'bottom');       % applied load

% Define the known parameters
sys = System(mesh);
sys.add_constant('P',10^3,'E',10^7,'a',10^-2);
sys.add_function('A',@(sys,elem,t) element_area(sys,elem));

% Define stiffness matrix with direct method
sys.add_matrix('K', 'A*E/L*Ke');

% Build the force vector
nat = mesh.get_dof('Boundary',3,'Component','x');
sys.add_vector('f',sys.get('P'),'dof',nat);

% Create solver, apply essential boundary conditions, and solve
solver = solvers.LinearSolver(sys);
solver.add_essential_boundary({'id',1,'value', 0},{'id',2,'Component','y','value',0});
u = solver.solve();

% Plot results
mesh.plot(u,'-deform','Component','x','colorbar','x-direction disp. (m)');

function A = element_area(sys,elem)
%ELEMENT_AREA returns the element area based on the element id

a = sys.get('a');
if elem.id == 2; 
    A = 2*a; 
else
    A = a; 
end



