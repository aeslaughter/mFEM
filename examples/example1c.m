%EXAMPLE1C Reproduces Example 5.1 of Fish & Belytschko (2007) using 
% automatic assembly of the stiffness matrix and force vector.
% 
% Syntax:
%   example1c
%   example1c('PropertyName', PropertyValue)
%
% Description:
%   example1c solves the propblem with 2 Line2 elements
%
%   example1c('PropertyName', PropertyValue) allows the user to customize
%   the behavior of this example using the property pairs listed below.
%
% Example1b Property Descriptions
%
% N
%   scalar
%   The number of elements in the x and y directions, the default is 2.
%
% Element
%   {'Line2'} | 'Line3'
%   Specifies the type of element for the mesh
%
% See also EXAMPLE1B
function example1c(varargin)

% Setup
import mFEM.* mFEM.solvers.*

% Set the default options and apply the user defined options
opt.n = 2;
opt.element = 'Line2';
opt = gather_user_options(opt,varargin{:});
    
% Create Mesh
mesh = FEmesh('Element', opt.element);
mesh.grid(0,4,opt.n);
mesh.init();

% Label the boundaries
mesh.add_boundary(1,'left');    % T = 0 boundary (essential)    
mesh.add_boundary(2,'right');   % q = 20 boundary 

% Build Matrix and Vector Equations
sys = System(mesh);
sys.add_constant('k',2,'A',0.1,'b',5,'q_bar',5);
sys.add_matrix('K', 'B''*k*A*B');
sys.add_vector('f', 'N''*b');
sys.add_vector('f', '-q_bar*A*N''', 'Boundary', 2);

solver = LinearSolver(sys);
solver.add_essential_boundary('id',1,'value',0);
T = solver.solve();

% figure;
% mesh.plot(T, '-ShowNodes');
% xlabel('Position, x (m)');
% ylabel('Temperature, T (C)');
% return;

% Compute the temperature gradients for the elements
k = 1;
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
mesh.plot(T,'Axes',h,'-ShowNodes'); hold on;
plot(h,x0,Tex,'k-','LineWidth',1);
legend({'FEM','Exact'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temperature (\circC)','interpreter','tex');

% Create TG plot
h = subplot(2,1,2);
plot(h,x0,TGex,'k-',TGx,TG,'bo','LineWidth',1);
legend({'Exact','FEM'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temp. Gradient (\circC/m)','interpreter','tex');
