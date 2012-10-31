function test
import mFEM.*

mesh = FEmesh();
mesh.add_element('Truss2',[-1,1; 0,0]);

N = mesh.element(1).shape(0);

N'*N