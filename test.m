function test
profile on;
import mFEM.*

mesh = FEmesh('Quad4','Type','CG');
mesh.grid(0,1,0,1,100,100);

% mesh.add_boundary('left',1);
% 
% mesh.plot();

% x = mesh.map.node;
% y = 2*x.^2;
% mesh.plot(y,'ElementLabels',true,'NodeLabels',true);
