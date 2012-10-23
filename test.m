function test
profile on;
import mFEM.*

mesh = FEmesh('Quad4');
mesh.grid(0,1,0,1,5,5);

elem = mesh.element(1);
elem. hmax()

% e = 1;
% for i = 1:4;
%     elem.side(i)
% end

