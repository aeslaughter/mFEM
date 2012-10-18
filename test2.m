function test2
clear;
import mFEM.*;

mesh = FEmesh('Quad4');
mesh.grid(0,1,0,1,1,1);

elem = mesh.element(1);

for i = 1:elem.n_sides();
    side = elem.build_side(i);
    n = side.get_normal()
end


