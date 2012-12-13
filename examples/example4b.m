% Example 9.2 of Fish & Belytschko (2007)
function example4b

% Import the mFEM library
import mFEM.*;
  
% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Element','Quad4','Space','vector');
mesh.add_element([0,1; 0,0; 2,0.5; 2,1]);
mesh.init();

% Label the boundaries
mesh.add_boundary(1, 'left');    % essential boundaries
mesh.add_boundary(2, 'top');     % distributed load (t = -20)

% Create system and add necessary components
sys = System(mesh);
sys.add_constant('E', 3e7, 'v', 0.3, 't_top', [0;-20]);
sys.add_constant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
sys.add_matrix('K', 'B''*D*B');
sys.add_vector('f', 'N''*t_top', 'Boundary', 2);

% Assemble the matrix and vector
K = sys.assemble('K'); 
f = sys.assemble('f');

% Assemble and solve
solver = solvers.LinearSolver(sys);
solver.add_essential_boundary('id',1,'value',0);
u = solver.solve()

% Compute the stress and strain at the Gauss points
% Loop through the elements
for e = 1:mesh.n_elements; % (include for illustration, but not needed)
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = u(elem.get_dof());
    
    % Compute the stress and strain at the Gauss points
    k = 1;
    qp = elem.quad.rules();
    for i = 1:length(qp);
        for j = 1:length(qp);
            strain(:,k) = elem.shape_deriv(qp(i),qp(j))*d;
            stress(:,k) = sys.get('D')*strain(:,k);
            k = k + 1;
        end
    end
end    

% Display the stress and strain vectors
strain, stress
