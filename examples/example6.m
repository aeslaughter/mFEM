% MAE4700/5700 HW8, Prob. 3
function example6

% Set known values
a = 5;
b = 10;
E = 10000;
u0 = -0.01;

% Import the mFEM library
import mFEM.*;
  
% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Space','vector','Element','Tri6');
mesh.grid(0, pi/2, a, b, 10, 5,'-pol2cart'); % x = theta; y = r
mesh.init();

% Label the boundaries
mesh.add_boundary(1, 'bottom');                     % y = 0 boundary (traction)
mesh.add_boundary(2, 'x < 0.001');                  % x = 0 boundary (rollers)
mesh.add_boundary(3,{'x < 0.001', ['y==',num2str(a)]});  % x = 0 and y = a (pin)

% Create system and add matrix components
sys = System(mesh);
sys.add_constant('E', E, 'v', 0.25);
sys.add_constant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
sys.add_matrix('K', 'B''*D*B');

% Assemble the matrix and vector
K = sys.assemble('K');

% Define dofs indices
ess(:,1) = mesh.get_dof('Boundary', 1, 'Component', 'x'); % traction
ess(:,2) = mesh.get_dof('Boundary', 2, 'Component', 'x'); % rollers
ess(:,3) = mesh.get_dof('Boundary', 3);                   % pin

% Define known displacements
u = zeros(mesh.n_dof,1);         % initialize the displacement vector
u(ess(:,1)) = u0;                % traction displacements
u(ess(:,2)) = 0;                 % zero displacements at rollers
u(ess(:,3)) = 0;                 % zero displacement at pin

% Combine essential boundaries for solution
ess = any(ess,2);

% Solve for the unknown displacements
ticID = tmessage('Solution step...');
u(~ess) = K(~ess,~ess)\(-K(ess,~ess)'*u(ess));
tmessage(ticID);

% Display the result
mesh.plot(u,'-new','component',1,'-deform','scale',100,'colorbar','u_x-displacement');
mesh.plot(u,'-new','component',2,'-deform','scale',100,'colorbar','u_y-displacement');
