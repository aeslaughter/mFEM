function test

import mFEM.*

mesh = FEmesh('Quad4','Type','DG');
mesh.grid(0,1,0,1,3,3);

e = 5;

elem = mesh.element(e);
for s = 1:elem.n_sides;
    elem.side(s).neighbor.id  
   elem.side(s) 
end

% mesh.plot();

