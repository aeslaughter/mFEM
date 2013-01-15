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
N = 16;
mesh.grid(0,1,0,1,N,N);
% mesh.grid(-1,0,0,1,N,N);
% mesh.grid(0,1,0,1,N,N);
% mesh.grid(-1,0,-1,0,N,N);
mesh.init();
% mesh.plot();

mesh.add_boundary(1);

sys = System(mesh);
sys.add_constant('p',-6);
sys.add_matrix('K','B''*B');
sys.add_vector('f','N''*p');

u_e = @(x) 1 + x(:,1).^2 + 2*x(:,2).^2;
solver = LinearSolver(sys);
solver.add_essential_boundary('id',1,'value',u_e);
u = solver.solve;
mesh.plot(u,'-new');

x = mesh.get_nodes();
ue = u_e(x);
mesh.plot(ue-u,'-new');






