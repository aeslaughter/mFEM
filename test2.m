function test2
clear all;


elem = Quad4(1,[0,1; 0,0; 2,0.5; 2,1]);
elem.side_detJ(4,0)


% addpath('..');
% mesh = FEmesh('Quad4');
% % mesh.gen2D(0,1,0,1,3,3);
% mesh.add_element([0,1; 0,0; 2,0.5; 2,1]);
% mesh.initialize();
% % mesh.plot();
% 
% mesh.side_shape


% mesh.add_boundary_id('top', 1);
% % mesh.add_boundary_id('right', 2);
% % mesh.add_boundary_id(3);
% 
% e = 1; s = 4;
% mesh.element(e).side(s)
% mesh.element(e).side(s).global_dof
% mesh.element(e).side(s).neighbor


end


