function test2



mesh = Mesh('Quad4','DG');
mesh.gen2D(0,1,0,1,2,1);
%mesh.map;

% elem = Quad4(1,x,y,'vector');
% 
% elem.shape(0,0);
% elem.shape_deriv(0,0)


% 
% elem.side_shape(1,0.5)
% 
% return
% 
% 
% N = @(xi,eta) elem.shape(xi,eta);
% J = @(xi,eta) elem.detJ(xi,eta);
% 
% N(-1,-1)
% N(1,-1)
% N(1,1)
% N(-1,1)
% 
% J(-1,1)
% J(1,-1)
% J(1,1)
% J(-1,1)
% 
% 
% 
% qp = [-1/sqrt(3), 1/sqrt(3)];
% 
% N(qp(1),qp(1))
% 
% N(qp(1),qp(1))' * N(qp(1),qp(1));