function example8

import mFEM.*

d_0 = @(x,y) (x-0.5).^2 + (y-0.75).^2 - 0.15^2;

mesh = FEmesh('Quad4','scalar','DG');
mesh.grid(0,1,0,1,1,1);

d = d_0(mesh.get_nodes(1), mesh.get_nodes(2));

% mesh.plot(d);

    function v = velocity(x,y,t)
       v(1,1) = cos(pi*t/8).*sin(2*pi*y).*sin(pi*x).^2; 
       v(2,1) = -cos(pi*t/8).*sin(2*pi*x).*sin(pi*y).^2; 
    end

q_elem = Gauss(2,'Quad');
[qp, W] = q_elem.rules();

M = sparse(mesh.n_dof, mesh.n_dof);
K = sparse(mesh.n_dof, mesh.n_dof);

t = 0;


for e = 1:mesh.n_elements;
    
   elem = mesh.element(e);
   
   Me = zeros(elem.n_dof);
   Ke = zeros(elem.n_dof);
   
   
   N = @(xi,eta) elem.shape(xi,eta);
   B = @(xi,eta) elem.shape_deriv(xi,eta);
   detJ = @(xi,eta) elem.detJ(xi,eta);
   
   for i = 1:length(qp);
       for j = 1:length(qp);
          
           
           [x,y] = elem.get_position(qp(i),qp(j));
           v = velocity(x,y,t);
           
           
           
           Me = Me + N(qp(i),qp(j))'*N(qp(i),qp(j))*detJ(qp(i),qp(j));
           Ke = Ke + (B(qp(i),qp(j))'*v)*N(qp(i),qp(j))*detJ(qp(i),qp(j));
           
   
       end
   end
   
   side = elem.build_side(1);
   GN = elem.local_grad_basis(0,-1);
   
   a = [GN(1,:)*elem.nodes,0]
   b = [GN(2,:)*elem.nodes,0]
   cross(a,b)
   
   GN*elem.nodes

%     Ke
    
    
end







end






