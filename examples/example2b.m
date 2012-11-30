% Example 8.1 of Fish & Belytschko (2007).
%
% Syntax
%   example2b
%
% Description
%   example2b solves a simple two element heat conduction problem.
function example2b
   
% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Element','Tri3');
mesh.add_element([0,0; 2,0.5; 0,1]);
mesh.add_element([2,0.5; 2,1; 0,1]);
mesh.init();

% Label the boundaries
mesh.add_boundary(1, 'top');     % q = 20 boundary
mesh.add_boundary(2, 'right');   % q = 0 boundary
mesh.add_boundary(3);            % essential boundaries (all others)

% Create the System
sys = System(mesh);
sys.add_constant('k', 5*eye(2), 'b', 6, 'q_top', 20);

% Create matrices
sys.add_matrix('K', 'B''*k*B');
sys.add_vector('f', 'N''*b');
sys.add_vector('f', 'N''*-q_top', 'Boundary', 1);

% Assemble
K = sys.assemble('K'); full(K)
f = sys.assemble('f');

% Define dof indices for the essential dofs and non-essential dofs
ess = mesh.get_dof('Boundary',3); % 4
non = ~ess;

% Solve for the temperatures
T = zeros(size(f));         % initialize the temperature vector
T(ess) = 0;                 % apply essential boundary condtions
T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

% Solve for the reaction fluxes
r = K*T - f;

% Display the results
T,r

% Compute the flux values for each element
% Loop through the elements
for e = 1:mesh.n_elements; 
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = T(elem.get_dof());
    
    % Compute the flux at the Gauss points
    q(:,e) = -sys.get('k')*elem.shape_deriv()*d;

end    

% Display the flux vectors
q