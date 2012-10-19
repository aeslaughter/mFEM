function Tri6_test

import mFEM.*

elem = Tri6(1, [0,0; 6,2; 4,4], 'Space', 'vector');

v = 1/3;
E = 288;
D = E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]; % constitutive matrix

q_elem = Gauss(2,'tri');
[qp, W] = q_elem.rules();

B = @(xi1, xi2) elem.shape_deriv(xi1,xi2);

Ke = zeros(elem.n_dof);
for i = 1:length(qp);
    Ke = Ke + W(i)*B(qp(i,1),qp(i,2))'*D*B(qp(i,1),qp(i,2))*elem.detJ(qp(i,1),qp(i,2));
end

Ke