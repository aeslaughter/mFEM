% MAE4700/5700 HW8, Prob. 3
function example8b

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
mesh.addBoundary(1, 'bottom');                     % y = 0 boundary (traction)
mesh.addBoundary(2, 'x < 0.001');                  % x = 0 boundary (rollers)
mesh.addBoundary(3,{'x < 0.001', ['y==',num2str(a)]});  % x = 0 and y = a (pin)

% Create system and add matrix components
sys = System(mesh);
sys.addConstant('E', E, 'v', 0.25);
sys.addConstant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
sys.addMatrix('K', 'B''*D*B');

% Assemble and Solve
solver = solvers.LinearSolver(sys,'force',0);
solver.addEssential({'boundary', 1, 'Component', 'x', 'value', u0},...
    {'boundary', 2, 'Component', 'x', 'value', 0},{'boundary', 3, 'value', 0});
u = solver.solve();

% Display the result
mesh.plot(u,'-new','component',1,'-deform','scale',100,'colorbar','u_x-displacement','Patch',{'EdgeColor','k'});
mesh.plot(u,'-new','component',2,'-deform','scale',100,'colorbar','u_y-displacement','Patch',{'EdgeColor','k'});
