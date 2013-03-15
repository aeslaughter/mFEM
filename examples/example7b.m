% MAE4700/5700 HW8, Ex. Prob.
function varargout = example7b(varargin)

% Gather inputs
opt.debug = false;
opt = gatherUserOptions(opt, varargin{:});


% Create a FEmesh object, add the single element, and initialize it
mesh = mFEM.Mesh('Space','vector');
L = 10;
c = 0.5;
mesh.grid('Quad4',0,L,-c,c,20,15);

% Label the boundaries
mesh.addBoundary(1, 'right'); 
mesh.addBoundary(2, 'left');
mesh.update();

% Create system and add matrix components
sys = mFEM.System(mesh);
sys.addConstant('c', c, 'l', L, 'E', 1e7, 'v', 0.3, 'P', [0;100]);
sys.addConstant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
sys.addMatrix('K', 'B''*D*B');
sys.addVector('f', 'N''*P', 'Tag', 2);

% Assemble the matrix and vector
solver = mFEM.solvers.LinearSolver(sys);
solver.addEssential('tag',1,'value',0);
u = solver.solve();

% Debug: no plot, return solution
if opt.debug;
    varargout = {u, exactSoln(sys)};
    return;
end

% Display the displacement results
subplot(2,1,1);
mesh.plot(u,'-deform','Patch',{'EdgeColor','k'},...
    'Colorbar','y-disp. (m)','Component', 2);
title('FEM Solution');
xlabel('x (m)'); 
set(gca,'Ytick',[]);
xlim([-0.1,10]);

% Display the exact solution
subplot(2,1,2);
u_exact = exactSoln(sys);
mesh.plot(u_exact,'-deform','Colorbar','y-disp. (m)','Component',2);
title('Exact Solution');
xlabel('x'); 
xlim([-0.1,10]);

function u = exactSoln(sys)
%EXACT_SOLN Specifies the analytical solution displacements

% Extract constants from System object
c = sys.get('c');  % 0.5;
E = sys.get('E');  % 1.0e7;
nu = sys.get('v'); % 0.3;
P = sys.get('P'); P = P(2); %100.0;
l = sys.get('l'); %10;

% Compute intermediate terms
I = 2.0/3*c^3;
G = E/(2*(1+nu));

% Define the exact solution function
u_x = @(x) -P*x(1)^2*x(2)/(2*E*I) - nu*P*x(2)^3/(6*E*I) + P*x(2)^3/(6*I*G) + (P*l^2/(2*E*I) - P*c^2/(2*I*G))*x(2);
u_y = @(x) nu*P*x(1)*x(2)^2/(2*E*I) + P*x(1)^3/(6*E*I) - P*l^2*x(1)/(2*E*I) + P *l^3/(3*E*I);

% Compute the exact solution
nodes = sys.mesh.getNodes();
x = nodes.getCoord();
u = zeros(2*length(x),1);
k = 1;
for i = 1:length(x);
    u(k) = u_x(x(i,:));
    u(k+1) = u_y(x(i,:));
    k = k + 2;
end
