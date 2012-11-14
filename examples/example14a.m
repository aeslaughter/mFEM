%EXAMPLE14a Beam element

function example14a

import mFEM.*

mesh = FEmesh('Element','Beam');
mesh.add_element([0;8]);
mesh.add_element([8;12]);
mesh.init();

mesh.add_boundary(1,'right');
mesh.add_boundary(2,'x==8');
mesh.add_boundary(3,'left');

P1 = -10;
P2 = 5;
EI = 10^4;
c = [-20,20]; % presecribed force and moment on boundary

K = Matrix(mesh);
f = zeros(mesh.n_dof,1);
for e = 1:mesh.n_elements;
    
    elem = mesh.element(e);
    
    [qp,w] = elem.quad.rules();
    B = @(xi) elem.shape_deriv(xi);
    N = @(xi) elem.shape(xi);
    
    Ke = zeros(elem.n_dof);    
    fe = zeros(elem.n_dof,1);
    for i = 1:length(qp);
        b = body_force(elem);
        Ke = Ke + w(i)*EI*B(qp(i))'*B(qp(i))*elem.detJ(qp(i));
        fe = fe + w(i)*N(qp(i))'*b*elem.detJ(qp(i));
    end

    if e == 1;
        fe = fe + N(0)'*P1;
    end
    
    if elem.boundary_id == 1;
        for s = 1:length(elem.side);
            if elem.side(s).boundary_id == 1;
                dof = elem.get_dof('Side',s,'-local');
                fe(dof) = c;
            end
        end
    end

    dof = elem.get_dof();
    K.add_matrix(Ke, dof);
    f(dof) = f(dof) + fe;
end
K = K.init();

dofP2 = mesh.get_dof('Boundary',2,'Component',1);
f(dofP2) = f(dofP2) + P2;

ess = mesh.get_dof('Boundary',3);
u = zeros(size(f));
u(ess) = 0;
u(~ess) = K(~ess,~ess)\f(~ess);

full(K)




function p = body_force(elem)
    x = elem.nodes(2);
    p = 0;
    if x <= 8;
        p = -1;
    end