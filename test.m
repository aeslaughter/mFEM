function test

clear;
import mFEM.*

mesh = FEmesh('Element','Quad4');
mesh.grid(0,1,0,1,3,3);
mesh.init();
% mesh.plot()

e = mesh.get_elements('contains',[0.4,0.4]);