function example8
    clear;
    import mFEM.*

    d_0 = @(x,y) (x-0.5).^2 + (y-0.75).^2 - 0.15^2;

%     if exist('example8.mat','file');
%         load('example8.mat');
%     else
        mesh = FEmesh('Type','DG');
        mesh.grid('Quad4',0,1,0,1,1,1);
        mesh.init();
        save('example8.mat','mesh');
%    end


    M = mass_matrix(mesh);  % inverse of mass matrix
    
    K = stiffness_matrix(mesh,0);
    full(K)
    
    
        return;
        
        
      % Initialize
    y0 = d_0(mesh.map.node(:,1), mesh.map.node(:,2));
    mesh.plot(y0);
  

    dt = 0.001;


    for t = 0:dt:2*dt%;

        K = stiffness_matrix(mesh,t);
        size(K)
        size(y0)

        y1 = y0 + dt*M*K*y0;
        y2 = 3/4*y0 + 1/4*y1 + 1/4*dt*M*K*y1;
        y0 = 1/3*y0 + 2/3*y2 + 2/3*dt*M*K*y2;

        mesh.plot(y0);



    end
    
end

function K = stiffness_matrix(mesh, t)

%
disp('Building stiffness matrix...');
tic;

import mFEM.*
q_elem = Gauss(2);
[qp, W] = q_elem.rules();

q_side = Gauss(1);
[qps, Ws] = q_side.rules();

K = Matrix(mesh);

n_ref = [1,0];

for e = 1:mesh.n_elements;
    
   elem = mesh.element(e);
   
   Ke = zeros(elem.n_dof);
   
   N = @(xi,eta) elem.shape(xi,eta);
   B = @(xi,eta) elem.shape_deriv(xi,eta);
   
   for i = 1:length(qp);
       for j = 1:length(qp);
           [x,y] = elem.get_position(qp(i),qp(j));
           v = velocity(x,y,t);
      
           Ke = Ke + W(i)*W(j)*(B(qp(i),qp(j))'*v)*N(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
       end
   end

   De = elem.get_dof;
   K.add_matrix(Ke, De);

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

            Np = @(beta) side_e.shape(beta);
            Nm = @(beta) side_n.shape(beta);

            de = elem.side(s_e).dof;
            dn = neighbor.side(s_n).dof;
            
            for i = 1:length(qps);
                [x,y] = side_e.get_position(qps(i));
                v = velocity(x,y,t);
                n = side_e.get_normal(qps(i));

                if dot(v,n) >= 0;
                    Ken(de,dn) = Ken(de,dn) - dot(n,n_ref)*Ws(i)*Np(qps(i))'*Nm(qps(i))*dot(v,n)*side_e.detJ();
                    Knn(dn,dn) = Knn(dn,dn) - dot(n,n_ref)*Ws(i)*Nm(qps(i))'*Nm(qps(i))*dot(v,n)*side_e.detJ();
                else
                    Kne(dn,de) = Kne(dn,de) - dot(n,n_ref)*Ws(i)*Nm(qps(i))'*Np(qps(i))*dot(v,n)*side_e.detJ();
                    Kee(de,de) = Kee(de,de) - dot(n,n_ref)*Ws(i)*Np(qps(i))'*Np(qps(i))*dot(v,n)*side_e.detJ();  
                end
            end

            De = elem.get_dof();
            Dn = neighbor.get_dof();
            
%             K.add_matrix(Kee, De, De);
%             K.add_matrix(Ken, De, Dn);
%             K.add_matrix(Kne, Dn, De);
%             K.add_matrix(Knn, Dn, Dn);
       end
   end 
end

% Compl
K = K.init();


disp(['   ...completed in ',num2str(toc),' sec.']);


end


function v = velocity(x,y,t)
   v(1,1) = cos(pi*t/8).*sin(2*pi*y).*sin(pi*x).^2; 
   v(2,1) = -cos(pi*t/8).*sin(2*pi*x).*sin(pi*y).^2; 
end

function M = mass_matrix(mesh)
import mFEM.*
q_elem = Gauss(2,'Quad');
[qp, W] = q_elem.rules();

M = Matrix(mesh.n_dof);

parfor e = 1:mesh.n_elements;
    
   elem = mesh.element(e);
   
   Me = zeros(elem.n_dof);

   
   N = @(xi,eta) elem.shape(xi,eta);
   
   for i = 1:length(qp);
       for j = 1:length(qp);
           Me = Me + W(i)*W(j)*N(qp(i),qp(j))'*N(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
       end
   end

   De = elem.get_dof;
   M.add_matrix(Me, De);
end

M = M.init();


end




