function test
profile on;
import mFEM.*

mesh = FEmesh('Quad4');
mesh.grid(0,1,0,1,50,50);

elem = mesh.element(1);
elem.size()

