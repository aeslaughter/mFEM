function test

import mFEM.*

mesh = FEmesh('Quad4');
mesh.grid(0,1,0,1,3,3);
% mesh.plot()

e = 5;
elem = mesh.element(e);

for i = 1:4;
    elem.side(i)
    elem.side(i).neighbor.id
end

