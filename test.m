function test
import mFEM.*

% q0 = Tri3(1,[0,1;1,0;1,1]);
% q1 = Tri3(2,[0,1;1,1;1,0]);
% any([q0,q1] == [q0])


mesh = FEmesh('Quad4','Type','DG');
mesh.grid(0,1,0,1,100,100);
% mesh.add_boundary('left',1);
% 
% mesh.plot();

% x = mesh.map.node;
% y = 2*x.^2;
% mesh.plot(y,'ElementLabels',true,'NodeLabels',true);
