function test

clear;
import mFEM.*

% mesh = FEmesh('Element','Hex8');
% mesh.grid(0,1,0,1,0,1,10,10,10);
% mesh.init();
% 
% % Initialize the temperatures
% T_exact = @(x,y,z,t) exp(-t)*sin(pi*x).*sin(pi*y).*sin(pi*z);
% nodes = mesh.get_nodes();
% T = T_exact(nodes(:,1),nodes(:,2),nodes(:,3),0);
% 
% mesh.plot(T);

mesh = FEmesh('Element','Quad4');
mesh.grid(0,1,0,1,10,10);
mesh.init();

mesh.add_subdomain(1,'x<0.5');

e = mesh.get_elements('Subdomain',1);


% Initialize the temperatures
% T_exact = @(x,y,t) exp(-t)*sin(pi*x).*sin(pi*y);
% nodes = mesh.get_nodes();
% T = T_exact(nodes(:,1),nodes(:,2),0);
% 
% mesh.plot(T);