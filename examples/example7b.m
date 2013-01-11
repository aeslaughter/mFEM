% MAE4700/5700 HW8, Ex. Prob.
function example7b

% Import the mFEM library
import mFEM.*;
  
% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Element','Quad4','Space','vector');
L = 10;
c = 0.5;
mesh.grid(0,L,-c,c,30,20);
mesh.init();

% Label the boundaries
mesh.add_boundary(1, 'right'); 
mesh.add_boundary(2, 'left');

% Create system and add matrix components
sys = System(mesh);
sys.add_constant('c', c, 'l', L, 'E', 1e7, 'v', 0.3, 'P', [0;100]);
sys.add_constant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
sys.add_matrix('K', 'B''*D*B');
sys.add_vector('f', 'N''*P', 'Boundary', 2);

% Assemble the matrix and vector
solver = solvers.LinearSolver(sys);
solver.add_essential_boundary('id',1,'value',0);
u = solver.solve();

% Display the displacement results
%subplot(2,1,1);
mesh.plot(u,'-deform','Patch',{'EdgeColor','k'},...
    'Colorbar','y-disp. (m)','Component', 2);
%title('FEM Solution');
xlabel('x (m)'); 
set(gca,'Ytick',[]);
return;
xlim([-0.1,10]);

% Display the exact solution
subplot(2,1,2);
u_exact = exact_soln(sys);
mesh.plot(u_exact,'-deform','Colorbar','y-disp. (m)','Component',2);
title('Exact Solution');
xlabel('x'); ylabel('y');
xlim([-0.1,10]);

function u = exact_soln(sys)
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
x = sys.mesh.get_nodes;
u = zeros(2*length(x));
k = 1;
for i = 1:length(x);
    u(k) = u_x(x(i,:));
    u(k+1) = u_y(x(i,:));
    k = k + 2;
end
