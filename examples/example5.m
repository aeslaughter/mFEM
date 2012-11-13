% MAE4700/5700 HW8, Prob. 2
function example5

% Import the mFEM library
import mFEM.*;
  
% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Space','vector');
mesh.grid('Quad4',0,2,0,1,1,1);
mesh.init();

% Label the boundaries
mesh.add_boundary(1, 'bottom','left');
mesh.add_boundary(2, 'right','top');

% Create system and add matrix components
sys = System(mesh);
sys.add_constant('E', 3e11, 'v', 0.3, 't', [1000;1000]);
sys.add_constant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
sys.add_matrix('K', 'B''*D*B');

% Add force components
sys.add_vector('f_1', 'N''*t', 'Boundary', 1);
sys.add_vector('f_2', 'N''*-t', 'Boundary', 2);

% Assemble the matrix and vector
K = sys.assemble('K'); 
f = sys.assemble('f_1') + sys.assemble('f_2');

% Define dof indices for the essential dofs and non-essential dofs
ess = [1,2,4];
non = [3,5,6,7,8];

% Solve for the temperatures
u = zeros(size(f));         % initialize the displacement vector
u(ess) = 0;                 % apply essential boundary condtions
u(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

% Solve for the reaction fluxes
r = K*u - f;

% Display the displacement results
u
mesh.plot(u,'-Deform','-ShowNodes','Colorbar','Magnitude of Disp. (m)');
xlabel('x'); ylabel('y');

% Compute the stress and strain at the Gauss points
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
