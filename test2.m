function test2
clear classes;
import mFEM.*;


e1 = Tri3(1,[0,0; 2,2; 0,4]);
e1.detJ()
e1.shape_deriv()


% e1 = Linear2(1,[0; 5]);
% e1.detJ()
% e1.shape_deriv()
% e1.side_detJ(2)


%e2 = Quad4(1,[0,1; 0,0; 2,0.5; 1,0]);
% e2.side_detJ(1,0)
% e2.side_detJ(2,0)
% e2.side_detJ(3,0)
% e2.side_detJ(4,0)


