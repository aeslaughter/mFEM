function test

clear;
import mFEM.*

mesh = FEmesh('Element','Quad4');
mesh.grid(0,1,0,1,3,3);
mesh.init();
% mesh.plot()

sys = System(mesh);

K = mFEM.base.Kernel(sys,'N''N');

