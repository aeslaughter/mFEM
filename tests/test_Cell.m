function test_Cell

import mFEM.cells.* mFEM.elements.base.*
% n(1) = Node(0);
% n(2) = Node(1);
c = Line([0,1],2);

c.build()