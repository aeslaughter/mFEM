function test

import mFEM.*

mesh = FEmesh('Quad4');
mesh.grid(0,1,0,1,100,100);
% mesh.plot()

elem = mesh.element(1);

