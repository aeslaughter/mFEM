% Example 8.1 of Fish & Belytschko (2007).
%
% Syntax
%   example2b
%
% Description
%   example2b solves a simple two element heat conduction problem.
function T = example2b
   
% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.addElement('Tri3',[0,0; 2,0.5; 0,1]);
mesh.addElement('Tri3',[2,0.5; 2,1; 0,1]);
mesh.init();

% Label the boundaries
mesh.addBoundary(1, 'top');     % q = 20 boundary
mesh.addBoundary(2, 'right');   % q = 0 boundary
mesh.addBoundary(3);            % essential boundaries (all others)

% Create the System
sys = System(mesh);
sys.addConstant('k', 5*eye(2), 's', 6, 'q_top', 20);

% Create matrices
sys.addMatrix('K', 'B''*k*B');
sys.addVector('f', 'N''*s');
sys.addVector('f', 'N''*-q_top', 'Boundary', 1);

% Assemble and solve
solver = solvers.LinearSolver(sys);
solver.addEssential('Boundary', 3, 'Value', 0);
T = solver.solve()

% Compute the flux values for each element
% Loop through the elements
for e = 1:mesh.n_elements; 
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = T(elem.getDof());
    
    % Compute the flux at the Gauss points
    q(:,e) = -sys.get('k')*elem.shapeDeriv()*d;

end    

% Display the flux vectors
q