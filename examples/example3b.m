% Example 8.2 of Fish & Belytschko (2007).
%
% Syntax:
%   example3
%
% Description:
%   example3 solves a simple single element heat conduction problem.
function example3b
   
% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.add_element('Quad4',[0,1; 0,0; 2,0.5; 2,1]);
mesh.init();

% Label the boundaries
mesh.add_boundary(1, 'top');     % q = 20 boundary
mesh.add_boundary(2, 'right');   % q = 0 boundary
mesh.add_boundary(3);       % essential boundaries (all others)

% Create the System
sys = System(mesh);
sys.add_constant('k', 5*eye(2), 'b', 6, 'q_top', 20);

% Create matrices
sys.add_matrix('K', 'B''*k*B');
sys.add_vector('f_s', 'N''*b');
sys.add_vector('f_q', 'N''*-q_top',1);

% Assemble
K = sys.assemble('K');
f = sys.assemble('f_s') + sys.assemble('f_q');

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

% Loop through the elements
for e = 1:mesh.n_elements; % (include for illustration, but not needed)
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = T(elem.get_dof());
    
    % Compute the flux at the Gauss points
    qp = elem.quad.rules();
    k = 1;
    for i = 1:length(qp);
        for j = 1:length(qp);
            q(:,k) = -sys.get('k')*elem.shape_deriv(qp(i),qp(j))*d;
            k = k + 1;
        end
    end
end    

% Display the flux vectors
q