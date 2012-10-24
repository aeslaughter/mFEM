function test_system

%% Setup
import mFEM.*

%% Create Mesh
mesh = FEmesh('Linear2');
mesh.grid(0,4,2);

mesh.add_boundary('left', 1); % T = 0 boundary (essential)    
mesh.add_boundary('right', 2); % q = 5 boundary  

%% Build Matrix and Vector Equations
sys = System(mesh);
sys.add_constant('k',2,'A',0.1,'b',5,'q_bar',5);
sys.add_matrix('K', 'B''*k*A*B');
sys.add_vector('f_s', 'N''*b');
sys.add_vector('f_q', '-q_bar*A*N''', 2);

K = sys.assemble('K');
f = sys.assemble('f_s') + sys.assemble('f_q');

%% Solve for the Temp. on non-essential boundaries
% Define dof indices for the essential dofs and non-essential dofs
non = mesh.get_dof(1,'ne'); % 2,3
ess = mesh.get_dof(1);      % 1

% Solve for the temperatures
T = zeros(size(f));         % initialize the temperature vector
T(ess) = 0;                 % apply essential boundary condtions
T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

%% Post Processing
% Compute the temperature gradients for the elements
% Loop through the elements
for e = 1:mesh.n_elements;
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Get the Gauss points
    qp = elem.quad.rules();
    
    % Collect the local values of T
    d(:,1) = T(elem.get_dof());
    
    % Compute the temperature gradient at the gauss point, store the value
    % twice for each element for creating graph, TGx is the node locations
    % used for plotting
    TG(1:2,e) = elem.shape_deriv(qp(1))*d;
    TGx(1:2,e) = elem.nodes;

end    

% Create exact solutions
x0 = 0:0.1:4;
Tex = -12.5*x0.^2 + 97.5*x0;
TGex = -25*x0 + 97.5;

% Generate figure for T and TG solutions
figure('Color','w','Name','Example 1 Results');

% Create Temperature plot
h = subplot(2,1,1);
x = unique(mesh.map.node);
plot(h,x0,Tex,'k-',x,T,'b-o','LineWidth',1);
legend({'Exact','FEM'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temperature (\circC)','interpreter','tex');

% Create TG plot
h = subplot(2,1,2);
plot(h,x0,TGex,'k-',TGx,TG,'b-o','LineWidth',1);
legend({'Exact','FEM'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temp. Gradient (\circC/m)','interpreter','tex');
