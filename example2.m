% Example 8.2 of Fish & Belytschko (2007)

function example2
    % 
    
clear classes;
    
addpath('..');
mesh = FEmesh('Quad4','vector');
mesh.add_element([0,1; 0,0; 2,0.5; 2,1]);
mesh.initialize();

mesh.add_boundary_id('left', 1);
mesh.add_boundary_id('top', 2);
mesh.add_boundary_id(3);

qrule = Gauss(2);
[qp, W] = qrule.rules();

qface = Gauss(1);
[qp_face, W_face] = qface.rules();

E = 3e7;
v = 0.3;

D = E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2];
tbar = [0;-20];

for e = 1:mesh.n_elements;
    
    elem = mesh.element(e);
    
    B = @(xi,eta) elem.shape_deriv(xi,eta);
    N = @(xi,eta) elem.shape(xi,eta);
    
    K = zeros(elem.n_dof);
    f = zeros(elem.n_dof,1);
    for i = 1:length(qp);
        for j = 1:length(qp);
%            B(qp(i),qp(j))
%             f = f + W(j)*W(i)*s*N(qp(i),qp(j))'*elem.detJ(qp(i),qp(j));
            K = K + W(j)*W(i)*B(qp(i),qp(j))'*D*B(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
        end
    end
    
    if any(elem.boundary_id == 1);
        for s = 1:elem.n_sides;
            if elem.side(s).boundary_id == 1;
                side = elem.build_side(s);
                
                for i = 1:length(qp_face);
                    N(-1,qp_face(i))'
                   f = f + (W_face(i)*side.shape(qp_face(i))'*side.detJ(qp_face(i)))*tbar;
                end
                delete(side);    
            end
        end          
    end
end

f
% dof = mesh.get_dof(3,'ne');
% T(4,1) = K(dof,dof)\f(dof);
% 
% r = K*T - f 




