function example8
clear;
import mFEM.*

d_0 = @(x,y) (x-0.5).^2 + (y-0.75).^2 - 0.15^2;

mesh = FEmesh('Quad4','scalar','DG');
mesh.grid(0,1,0,1,3,3);
%mesh.plot()


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
return;

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


   for i = 1:1%elem.n_neigbors;
        neighbor = elem.neighbors(i);

        
        Kee = zeros(elem.n_dof);
        Ken = zeros(elem.n_dof, neighbor.n_dof);
        Kne = zeros(neighbor.n_dof, elem.n_dof);
        Knn = zeros(neighbor.n_dof);
       
        Nn = @(xi,eta) neighbor.shape(xi,eta);
        
        [se, sn] = elem.match_sides(neighbor);
        
        side_e = elem.build_side(se);
        side_n = neighbor.build_side(sn); 
        
        N_plus = @(beta) side_e.shape(beta);
        N_minus = @(beta) side_n.shape(beta);

   end
   
   
    
end







end






