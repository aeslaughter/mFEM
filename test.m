function test

import mFEM.*

elem = Tri6(1,[0,0; 1,0; 0,1]);

N = elem.shape(0.5,0.5);
B = elem.shape_deriv(0.5,0.5)
