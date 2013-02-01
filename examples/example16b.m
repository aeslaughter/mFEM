% Solves Possion equation
%
% Strong Form:
%   -\nabla^2u - p = 0
% 
% Solution:
%   p = -6, u = 1 + x^2 + 2*y^2
%
% Syntax:
%   example16b
%
% Description:
%   example16b solves 2D Possion equation
function example16b
   
% Import the mFEM library
import mFEM.* mFEM.solvers.*

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Element','Quad4','Type','CG');
N = 10;
mesh.grid(0,1,0,1,N,N);
mesh.grid(-1,0,0,1,N,N);
mesh.grid(-1,0,-1,0,N,N);
mesh.init();

mesh.addBoundary(1);

sys = System(mesh);
sys.addConstant('p',-6);
sys.addMatrix('K','B''*B');
sys.addVector('f','N''*p');

u_e = @(x) 1 + x(:,1).^2 + 2*x(:,2).^2;
solver = LinearSolver(sys);
solver.addEssential('boundary',1,'value',u_e);
u = solver.solve;
mesh.plot(u,'-new');

x = mesh.getNodes();
ue = u_e(x);
mesh.plot(u_e(x),'-new');
mesh.plot(ue-u,'-new');






