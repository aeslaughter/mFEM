function example8
clear;
import mFEM.*

d_0 = @(x,y) (x-0.5).^2 + (y-0.75).^2 - 0.15^2;

mesh = FEmesh('Quad4','Type','DG');
mesh.grid(0,1,0,1,2,2);
% mesh.plot()


d = d_0(mesh.get_nodes(1), mesh.get_nodes(2));

% mesh.plot(d);

    function v = velocity(x,y,t)
       v(1,1) = cos(pi*t/8).*sin(2*pi*y).*sin(pi*x).^2; 
       v(2,1) = -cos(pi*t/8).*sin(2*pi*x).*sin(pi*y).^2; 
    end

q_elem = Gauss(2,'Quad');
[qp, W] = q_elem.rules();

q_side = Gauss(1,'Quad');
[qp_side, W_side] = q_side.rules();

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

   De = elem.get_dof;
   K(De,De) = K(De,De) + Ke;

   for s_e = 1:elem.n_sides;
       if ~elem.side(s_e).on_boundary;

            neighbor = elem.side(s_e).neighbor;
            s_n = elem.side(s_e).neighbor_side;

            Kee = zeros(elem.n_dof);
            Ken = zeros(elem.n_dof, neighbor.n_dof);
            Kne = zeros(neighbor.n_dof, elem.n_dof);
            Knn = zeros(neighbor.n_dof);

            side_e = elem.build_side(s_e);
            side_n = neighbor.build_side(s_n); 

            N_plus = @(beta) side_e.shape(beta);
            N_minus = @(beta) side_n.shape(beta);

            de = elem.side(s_e).dof;
            dn = neighbor.side(s_n).dof;
            
            for i = 1:length(qp_side);
                [x,y] = side_e.get_position(qp_side(i));
                v = velocity(x,y,t);
                n = side_e.get_normal(qp_side(i));
                
                if dot(v,n) >= 0;
                    Ken(de,dn) = Ken(de,dn) + N_plus(qp_side(i))'*N_minus(qp_side(i))*dot(v,n)*side_e.detJ();
                    Knn(dn,dn) = Knn(dn,dn) + N_minus(qp_side(i))'*N_minus(qp_side(i))*dot(v,n)*side_e.detJ();
                else
                    Kne(dn,de) = Kne(dn,de) + N_minus(qp_side(i))'*N_plus(qp_side(i))*dot(v,n)*side_e.detJ();
                    Kee(de,de) = Kee(de,de) + N_plus(qp_side(i))'*N_plus(qp_side(i))*dot(v,n)*side_e.detJ();  
                end
                    
            end

            De = elem.get_dof();
            Dn = neighbor.get_dof();
            
            K(De,De) = K(De,De) + Kee;
            K(De,Dn) = K(De,Dn) + Ken;            
            K(Dn,De) = K(Dn,De) + Kne;
            K(Dn,Dn) = K(Dn,Dn) + Knn;
       end
   end 
end



end






