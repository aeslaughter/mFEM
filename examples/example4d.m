% MAE4700/5700 HW8, Prob. 3
function example4d

% Set known values
a = 5;
b = 10;
E = 10000;
u0 = -0.01;

% Import the mFEM library
import mFEM.*;
  
% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Space','vector');
mesh.grid('Quad4', 0, pi/2, a, b, 10, 5); % x = theta; y = r
mesh.init();

% Label the boundaries
mesh.add_boundary('left', 1); 
mesh.add_boundary('right', 2);
% mesh.add_boundary({'x==pi/2','y==a'},2);

% Create system and add matrix components
sys = System(mesh);
sys.add_constant('E', E, 'v', 0.25);
sys.add_constant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
sys.add_matrix('K', 'B''*D*B');

% Assemble the matrix and vector
K = sys.assemble('K');

% Define dof indices for the essential dofs 
% (fixed displacement in theta direction)
ess(:,1) = mesh.get_dof('Boundary', 1, 'Component', 'x');

% Locate theta = pi/2 and r = a location (this will be automatic)
idx = mesh.map.node(:,1) == pi/2 & mesh.map.node(:,2) == a;
dof = mesh.map.dof(idx);
dof = transform_dof(dof, mesh.n_dof_node);
ess(dof,1) = true;

% Define dofs indices for u_0 (displacement in r-direction)
ess(:,2) = mesh.get_dof('Boundary', 2, 'Component', 'y');

% Define known displacements
u = zeros(mesh.n_dof,1);         % initialize the displacement vector
u(ess(:,1)) = 0;                 % fixed displacements
u(ess(:,2)) = u0;                % traction displacements

% Combine essential boundaries for solution
ess = any(ess,2);

% Solve for the unknown displacements
u(~ess) = K(~ess,~ess)\(-K(ess,~ess)'*u(ess));
% u(~ess) = K(~ess,~ess)\(f(~ess) - K(ess,~ess)'*u(ess));

% Display the results
mesh.plot(u,'polar',true,'component',2);

% Display the displacement results
% subplot(2,1,1);
% mesh.plot(u,'component',1,'colorbar','theta-displacement');
% xlabel('theta'); ylabel('r');
% 
% subplot(2,1,2);
% mesh.plot(u,'component',2,'colorbar','r-displacement');
% xlabel('theta'); ylabel('r');


