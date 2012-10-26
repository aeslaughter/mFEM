function test
import mFEM.*

mesh = FEmesh();
% mesh.add_element('Tri6',[0,0; 1,0; 0,1]);
mesh.grid('Tri6',0,1,0,1,10,10);
mesh.init();
% mesh.plot();

[mesh.map.node,mesh.map.dof]
% node = unique(mesh.map.node,'rows','R2012a');
x = mesh.map.node(:,1);
y = mesh.map.node(:,2);
T = @(x,y,t) exp(-t)*sin(pi*x).*sin(pi*y);
TT = T(x,y,0)
mesh.plot(TT(mesh.get_dof()));



% mesh.plot();

% elem = mesh.element(1);
% idx = elem.node_plot_order
% elem.nodes
% elem.nodes(idx,:)

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
