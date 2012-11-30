function example11

% Import the finite element library
import mFEM.*

% Set the initial conditions for the signed distance function
d_0 = @(x,y) (x-0.5).^2 + (y-0.75).^2 - 0.15^2;

% Create the mesh
mesh = FEmesh('Type','DG');
mesh.grid('Quad4',0,1,0,1,10,10);
mesh.init();

% Extract the x,y coordinates from the FEmesh and build a grid version
x = mesh.map.node(:,1);
y = mesh.map.node(:,2);
[X,Y] = meshgrid(x,y);

% Compute the intial condition
y0 = d_0(mesh.map.node(:,1), mesh.map.node(:,2));

% Compute the inverse of the mass matrix (constant w/ time)
M = mass_matrix(mesh); 

% Graph the initial condition
mesh.plot(y0,'-new'); hold on;
Z = griddata(x,y,y0,X,Y);
[~,chandle] = contour(X,Y,Z,[0,0],'-w','LineWidth',2);
title('t=0');
drawnow;

% Compute terms for initial time step
h = mesh.element(1).hmin();
c = 0;
for i = 1:length(x);
    v = velocity(x(i),y(i),0);
    if norm(v) > c;
        c = norm(v); 
    end
end

% Perform temporal calculations
t = 0;
while t <= 1;

    dt = 0.5 * h/(c*(2*3+1));
    t = t + dt;
    
    str = sprintf('Time %g (dt = %g)',t,dt);
    disp(str);
    
    [K,c] = stiffness_matrix(mesh,t);
        
    y1 = y0 + dt*M*K*y0;
    y2 = 3/4*y0 + 1/4*y1 + 1/4*dt*M*K*y1;
    y0 = 1/3*y0 + 2/3*y2 + 2/3*dt*M*K*y2;

    mesh.plot(y0); hold on;
    Z = griddata(x,y,y0,X,Y);
    delete(chandle);
    [~,chandle] = contour(X,Y,Z,[0,0],'-w','LineWidth',2);
    title(['t = ',num2str(t)]);
    drawnow;
end
    
end

function [K,v_max] = stiffness_matrix(mesh, t)

v_max = 0;

ticID = tmessage('Assembling stiffness matrix...');

import mFEM.*

K = Matrix(mesh);

n_ref = [1,0];

for e = 1:mesh.n_elements;
    
   elem = mesh.element(e);
      [qp, W] = elem.quad.rules('-cell');

   Ke = zeros(elem.n_dof);
   
   N = @(xi,eta) elem.shape(xi,eta);
   B = @(xi,eta) elem.shape_deriv(xi,eta);
   
   for i = 1:length(qp);
       [x,y] = elem.get_position(qp{i,:});
       v = velocity(x,y,t);
       Ke = Ke + W(i)*(B(qp{i,:})'*v)*N(qp{i,:})*elem.detJ(qp{i,:});

       if norm(v) > v_max;
           v_max = norm(v);
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

            de = elem.get_dof('Side',s_e,'-local');
            dn = neighbor.get_dof('Side',s_n,'-local');
            
            [qp,W] = side_e.quad.rules('-cell');
            
            for i = 1:length(qp);
                [x,y] = side_e.get_position(qp{i,:});
                v = velocity(x,y,t);
                n = side_e.get_normal(qp{i,:});

                if dot(v,n) >= 0;
                    Ken(de,dn) = Ken(de,dn) - W(i)*Np(qp{i,:})'*Nm(qp{i,:})*dot(v,n)*side_e.detJ();
                    Knn(dn,dn) = Knn(dn,dn) - W(i)*Nm(qp{i,:})'*Nm(qp{i,:})*dot(v,n)*side_e.detJ();
                else
                    Kne(dn,de) = Kne(dn,de) - W(i)*Nm(qp{i,:})'*Np(qp{i,:})*dot(v,n)*side_e.detJ();
                    Kee(de,de) = Kee(de,de) - W(i)*Np(qp{i,:})'*Np(qp{i,:})*dot(v,n)*side_e.detJ();  
                end
            end

            De = elem.get_dof();
            Dn = neighbor.get_dof();
            
            K.add_matrix(Kee, De, De);
            K.add_matrix(Ken, De, Dn);
            K.add_matrix(Kne, Dn, De);
            K.add_matrix(Knn, Dn, Dn);
       end
   end 
end

% Compl
K = K.init();

tmessage(ticID);


end

function v = velocity(x,y,t)
   v(1,1) = cos(pi*t/8).*sin(2*pi*y).*sin(pi*x).^2; 
   v(2,1) = -cos(pi*t/8).*sin(2*pi*x).*sin(pi*y).^2; 
end

function M = mass_matrix(mesh)
import mFEM.*


M = Matrix(mesh.n_dof);

for e = 1:mesh.n_elements;
    
   elem = mesh.element(e);
   [qp, W] = elem.quad.rules('-cell');
   
   Me = zeros(elem.n_dof);

   N = @(xi,eta) elem.shape(xi,eta);
   
   for i = 1:length(qp);
        Me = Me + W(i)*N(qp{i,:})'*N(qp{i,:})*elem.detJ(qp{i,:});
   end

   De = elem.get_dof;
   M.add_matrix(inv(Me), De);
end

M = M.init();

end
