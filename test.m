function test
import mFEM.*

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.grid('Quad4', 0, pi/2, 5, 10, 10, 10,'-pol2cart'); % x = theta; y = r
mesh.init();
mesh.plot();



% mesh = FEmesh('Space','vector');
% mesh.grid('Quad4',0,1,0,1,2,2);
% mesh.init();
% 
% 
% mesh.add_boundary(1,{'x==1','y==1'})
% 

