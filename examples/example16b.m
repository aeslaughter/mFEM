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
mesh = FEmesh('Element','Quad4','Type','CG','-time');
N = 10;
mesh.grid(0,1,0,1,N,N);
mesh.grid(-1,0,0,1,N,N);
mesh.grid(-1,0,-1,0,N,N);
mesh.init();

mesh.addBoundary(1);

sys = System(mesh,'-time');
sys.addConstant('p',-6);
sys.addMatrix('K','B''*B');
sys.addVector('f','N''*p');

u_e = @(x) 1 + x(:,1).^2 + 2*x(:,2).^2;
%u_e = @(x) exact(x);
solver = LinearSolver(sys);
solver.addEssential('boundary',1,'value',u_e);
u = solver.solve;
mesh.plot(u,'-new');
title('FE Solution');

x = mesh.getNodes();
mesh.plot(u_e(x),'-new');
title('Exact Solution');

mesh.plot(u_e(x)-u,'-new');
title('Error');

function u = exact(x)
[theta,r] = cart2pol(x(:,1),x(:,2));
if theta< 0; theta = theta + 2*pi(); end
u = r.^(2/3).*sin(2/3*theta);




