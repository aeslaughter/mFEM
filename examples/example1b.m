%EXAMPLE1B Reproduces Example 5.1 of Fish & Belytschko (2007) using 
% automatic assembly of the stiffness matrix and force vector.
% 
% Syntax:
%   example1b
%   example1b('PropertyName', PropertyValue)
%
% Description:
%   example1b solves the propblem with 2 Line2 elements
%
%   example1b('PropertyName', PropertyValue) allows the user to customize
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
function T = example1b(varargin)

% Set the default options and apply the user defined options
opt.n = 2;
opt.element = 'Line2';
opt = gatherUserOptions(opt,varargin{:});
    
% Create Mesh
mesh = mFEM.Mesh();
mesh.grid(opt.element,0,4,opt.n);
mesh.init();

% Label the boundaries
mesh.addBoundary(1,'left');    % T = 0 boundary (essential)    
mesh.addBoundary(2,'right');   % q = 20 boundary 

% Build Matrix and Vector Equations
sys = mFEM.System(mesh);
sys.addConstant('k',2,'A',0.1,'b',5,'q_bar',5);
sys.addMatrix('K','B''*k*A*B');
sys.addVector('f','N''*b');
sys.addVector('f','-q_bar*A*N''','Tag',2);

% Assemble and solve
solver = mFEM.solvers.LinearSolver(sys);
solver.addEssential('Tag',1,'Value',0);
T = solver.solve();

% Compute the temperature gradients for the elements
k = 1;
elements = mesh.getElements();
for e = 1:length(elements);
    
    % Extract the current element from the mesh object
    elem = elements(e);
    
    % Get the Gauss points
    qp = elem.quad.rules();
    
    % Collect the local values of T
    d(:,1) = T(elem.getDof());
    
    % Compute the temperature gradient at the gauss point, store the value
    % twice for each element for creating graph, TGx is the node locations
    % used for plotting
    for i = 1:length(qp);
        TG(k) = elem.shapeDeriv(qp(i))*d;
        TGx(k) = elem.getPosition(qp(i));
        k = k + 1;
    end
end    

% Build exact solutions
x0 = 0:0.1:4;
Tex = -12.5*x0.^2 + 97.5*x0;
TGex = -25*x0 + 97.5;

% Create Temperature plot
h = subplot(2,1,1);
mesh.plot(T,'Axes',h,'-showNodes'); hold on;
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
