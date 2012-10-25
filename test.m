function test
import mFEM.*

mesh = FEmesh('Linear2');
mesh.add_element([0;2]);
mesh.initialize();

elem = mesh.element(1);

side = elem.build_side(1);
side.shape()


% q = Gauss(3,'tri')
% [qp,w] = q.rules('cell')
% 
% qp{1}(:)
% 


% q0 = Tri3(1,[0,1;1,0;1,1]);
% 
% nargin('mFem.Tri3')
% 
% % q1 = Tri3(2,[0,1;1,1;1,0]);
% % any([q0,q1] == [q0])
% 
% 
% mesh = FEmesh('Quad4','Type','DG');
% mesh.grid(0,1,0,1,100,100);



% mesh.add_boundary('left',1);
% 
% mesh.plot();

% x = mesh.map.node;
% y = 2*x.^2;
% mesh.plot(y,'ElementLabels',true,'NodeLabels',true);
