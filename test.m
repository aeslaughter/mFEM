function test

clear;
import mFEM.*

% elem = Line2(1,[0,0;1,1],'space','vector');
% 
% x = elem.get_position(1)

mesh = FEmesh('Element','Hex8');
mesh.grid(0,1,0,1,0,1,1,1,1);
mesh.init();

mesh.plot();


side = mesh.element(1).build_side(2)