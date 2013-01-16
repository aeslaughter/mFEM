% Solves Possion equation
%
% Strong Form:
%   -\nabla^2u - p = 0
% 
% Solution:
%   p = -6, u = 1 + x^2 + 2*y^2
%
% Syntax:
%   example17b
%
% Description:
%   example17b solves 2D Possion equation
function example17b
   
% Import the mFEM library
import mFEM.* mFEM.solvers.*


% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Element','Quad4','Type','CG');
N = 8;
mesh.grid(-1,0,0,1,N,N);
mesh.grid(-1,0,-1,0,N,N);
mesh.grid(0,1,0,1,N,N);
mesh.init();
% mesh.plot();

mesh.add_boundary(1);

% sys = System(mesh);
% sys.add_constant('p',0);
% sys.add_matrix('K','B''*B');
% sys.add_vector('f','N''*p');

%u_e = @(x) x(:,1).^(2/3).*sin(2/3*x(:,2));
%u_e = @(x) cos(x(:,1)) .* exp(x(:,2));
% solver = LinearSolver(sys);
% solver.add_essential_boundary('id',1,'value', @(x) exact(x));
% u = solver.solve;
% mesh.plot(u,'-new');

x = mesh.get_nodes();
ue = exact(x);
mesh.plot(ue,'-new');


function u = exact(x)
[theta,r] = cart2pol(x(:,1),x(:,2));
% if theta < 0; theta = theta + 2*pi(); end
u = r.^(2/3).*sin(2/3*theta);



