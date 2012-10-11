% Example 8.2 of Fish & Belytschko (2007) example1
function example1
   
% Clear all variables, including classes and add the source directory
clear classes;
import mFEM.*;

mesh = FEmesh('Quad4');
mesh.add_element([0,1; 0,0; 2,0.5; 2,1]);
mesh.initialize();

mesh.add_boundary_id('top', 1);
mesh.add_boundary_id('right', 2);
mesh.add_boundary_id(3);

qrule = Gauss(2);
[qp, W] = qrule.rules();

qface = Gauss(1);
[qp_face, W_face] = qface.rules();

k = 5*eye(2);
s = 6;
q = 20;



for e = 1:mesh.n_elements;
    
    elem = mesh.element(e);
    
    B = @(xi,eta) elem.shape_deriv(xi,eta);
    N = @(xi,eta) elem.shape(xi,eta);
    N_side = @(id, beta) elem.side_shape(id, beta);
    
    K = zeros(elem.n_dof);
    f = zeros(elem.n_dof,1);
    for i = 1:length(qp);
        for j = 1:length(qp);
            f = f + W(j)*W(i)*s*N(qp(i),qp(j))'*elem.detJ(qp(i),qp(j));
            K = K + W(j)*W(i)*B(qp(i),qp(j))'*k*B(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
        end
    end
    
    
    for s = 1:elem.n_sides;
        if elem.side(s).boundary_id == 1;
            for i = 1:length(qp_face);
               f = f + -q*W_face(i)*N_side(s,qp_face(i))'*elem.side_detJ(s,qp_face(i));
               elem.side_detJ(s,qp_face(i))
               
            end
        end
    end          
end

f
dof = mesh.get_dof(3,'ne');
T(4,1) = K(dof,dof)\f(dof);

r = K*T - f 




