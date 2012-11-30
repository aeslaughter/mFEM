function test

clear;
import mFEM.*

mesh = FEmesh('Element','Beam');
mesh.grid(0,1,1);
mesh.init();

elem = mesh.element(1);
Ke = elem.stiffness()

[qp,w] = elem.quad.rules();
Ke = zeros(elem.n_dof);
for i = 1:length(qp);
    Ke = Ke + w(i)*elem.shape_deriv(qp(i))'*elem.shape_deriv(qp(i))*elem.detJ(qp(i));
end
Ke