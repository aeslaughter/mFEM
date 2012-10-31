% Example1c
% Reproduces Example 5.1 of Fish & Belytschko (2007) using automatic 
% assembly of the stiffness matrix and force vector and Line3 element.
% 
% Syntax:
%   example1c
%   example1c(n)
%
% Description
%   example1c runs example exactly as done in textbook, with two elements
%   example1c(n) runs the example with n number of elements
%
% See also EXAMPLE1B
function example1c(varargin)

% Setup
import mFEM.*

% Parse optional input
nel = 2;
if nargin == 1 && isnumeric(varargin{1}); 
    nel = varargin{1};
end
    
% Create Mesh
mesh = FEmesh();
mesh.grid('Line3',0,4,nel);
mesh.init();

mesh.add_boundary('left', 1); % T = 0 boundary (essential)    
mesh.add_boundary('right', 2); % q = 5 boundary  

% Build Matrix and Vector Equations
sys = System(mesh);
sys.add_constant('k',2,'A',0.1,'b',5,'q_bar',5);
sys.add_matrix('K', 'B''*k*A*B');
sys.add_vector('f_s', 'N''*b');
sys.add_vector('f_q', '-q_bar*A*N''', 2);

K = sys.assemble('K');
f = sys.assemble('f_s') + sys.assemble('f_q');

% Solve for the Temp. on non-essential boundaries
% Define dof indices for the essential dofs and non-essential dofs
non = mesh.get_dof(1,'ne');
ess = mesh.get_dof(1);

% Solve for the temperatures
T = zeros(size(f));         % initialize the temperature vector
T(ess) = 0;                 % apply essential boundary condtions
T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

% Compute the temperature gradients for the elements
k = 1;
for e = 1:mesh.n_elements; % loop through the elements
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Get the Gauss points
    qp = elem.quad.rules();
    
    % Collect the local values of T
    d(:,1) = T(elem.get_dof());
    
    % Compute the temperature gradient at the gauss point, store the value
    % twice for each element for creating graph, TGx is the node locations
    % used for plotting
    for i = 1:length(qp);
        TG(k) = elem.shape_deriv(qp(i))*d;
        TGx(k) = elem.get_position(qp(i));
        k = k + 1;
    end
end    

% Build exact solutions
x0 = 0:0.1:4;
Tex = -12.5*x0.^2 + 97.5*x0;
TGex = -25*x0 + 97.5;

% Create Temperature plot
h = subplot(2,1,1);
mesh.plot(T,'AxesHandle',h,'ShowNodes','on'); hold on;
plot(h,x0,Tex,'k-','LineWidth',1);
legend({'FEM','Exact'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temperature (\circC)','interpreter','tex');

% Create TG plot
h = subplot(2,1,2);
plot(h,x0,TGex,'k-',TGx,TG,'b-o','LineWidth',1);
legend({'Exact','FEM'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temp. Gradient (\circC/m)','interpreter','tex');
