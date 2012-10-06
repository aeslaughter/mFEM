function test

% Velocity field
    function u = velocity(x,y,t)
        u(1,1) = cos(pi*t/8)*sin(2*pi*y)*sin(pi*x)^2;
        u(2,1) = cos(pi*t/8)*sin(2*pi*x)*sin(pi*y)^2;
    end

% Phi initial
    function d = phi_initial(x,y)
        d = (x - 0.5)^2 + (y - 0.75)^2 - 0.15^2;
    end

% x = (0:0.05:1)';
% [X,Y] = meshgrid(x,x);
% figure('NextPlot','replace');
% %for t = 0:0.1:8;
%     for i = 1:length(x); 
%         for j = 1:length(x);
%            %u(i,j,:) = velocity(X(i,j),Y(i,j),t);
%            u(i,j) = phi_initial(X(i,j),Y(i,j));
%         end;
%     end;
%     pause(1);
%     contour(X,Y,u);
% %end
% 
% return;


mesh = Mesh('Quad4');
mesh.gen2D(0,1,0,1,2,2);

ndof = mesh.n_elements*4;
M = sparse(ndof, ndof);
K = sparse(ndof, ndof);

D = 5*eye(2);
qrule = Gauss(2);


t = 0;


for e = 1:mesh.n_elements;
    
    elem = mesh.element(e);
    
    N = @(xi,eta) elem.shape(xi,eta);
   % B = @(xi,eta) elem.shape_deriv(xi,eta);

     Me = zeros(elem.n_shape);
    %Ke = zeros(elem.n_shape);

    [qp, W] = qrule.rules();

    for i = 1:length(qp);
        for j = 1:length(qp);
           % v = velocity(qp(i), qp(j), t)
           % dot(v*N(qp(i),qp(j)), B(qp(i),qp(j)))

           % dot(v,B(qp(i),qp(j)))
            
            Me = Me + W(i)*N(qp(i),qp(j))'*N(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
            
            
            
           % Ke = Ke + W(i)*N(qp(i),qp(j))'*B(qp(i),qp(j)))*elem.detJ(qp(i),qp(j));
        end
    end
    
    dof = elem.global_dof;
    
    M(dof, dof) = Me;
    %K(dof, dof) = Ke;

end

%K

end
