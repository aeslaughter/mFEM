function test
profile on;
import mFEM.*

mesh = FEmesh('Quad4');
mesh.grid(0,1,0,1,100,100);
% mesh.plot()
% profile viewer


% 
% 
% e = 5;
% elem = mesh.element(e);
% for i = 1:4;
%     elem.side(i).neighbor.id
%     elem.side(i)
% end

