function example9
   
% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.add_element('Truss2',[-1,1; 0,0]);
mesh.add_element('Truss2',[0,1; 0,0]);
mesh.add_element('Truss2',[1,1; 0,0]);
mesh.init();

% Label the boundaries
mesh.add_boundary('top', 1);     % essential boundaries (fixed)

% Create the System
sys = System(mesh);
sys.add_constant('E', 10e7, 'A1', 10e-2, 'A2', '2*A1', 'A3', 'A1');

% Create matrices
sys.add_matrix('K_1', 'B''*A1*E*B');

K_1 = sys.assemble('K_1');
full(K_1)

% sys.add_vector('f_s', 'N''*b');
% sys.add_vector('f_q', 'N''*-q_top',1);
% 
% % Assemble
% K = sys.assemble('K');
% f = sys.assemble('f_s') + sys.assemble('f_q');
% 
% % Define dof indices for the essential dofs and non-essential dofs
% non = mesh.get_dof(3,'ne'); % 4
% ess = mesh.get_dof(3);      % 1,2,3
% 
% % Solve for the temperatures
% T = zeros(size(f));         % initialize the temperature vector
% T(ess) = 0;                 % apply essential boundary condtions
% T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries
% 
% % Solve for the reaction fluxes
% r = K*T - f;
% 
% % Display the results
% T,r