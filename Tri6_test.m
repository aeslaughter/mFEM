% A test function for the Tri6 Element
%
% Computes the element stiffness matrix for two test cases. The element
% equations and the test cases are from:
% http://mmc.geofisica.unam.mx/Bibliografia/Matematicas/EDP/MetodosNumericos/FEM/IFEM.Ch24.pdf

function Tri6_test

% Importh the FEM library
import mFEM.*

% Test 1 Element
elem = Tri6(1,[0,0; 6,2; 4,4], 'Space', 'vector');

% Define constitutive properties
v = 1/3;
E = 288;
D = E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]; % constitutive matrix

% Get the Gauss quadrature rules
q_elem = Gauss(2,'tri');
[qp, W] = q_elem.rules();

% Short-hand for shape function derivatives
B = @(xi1, xi2) elem.shape_deriv(xi1,xi2);

% Compute element stiffness matrix (compare to Eq. 24.35)
Ke = zeros(elem.n_dof);
for i = 1:length(qp);
    Ke = Ke + W(i)*B(qp(i,1),qp(i,2))'*D*B(qp(i,1),qp(i,2))*elem.detJ(qp(i,1),qp(i,2));
end
Ke

% Test2 Element
elem = Tri6(2,[-1/2,0; 1/2,0; 0, sqrt(3)/2; 0, -1/(2*sqrt(3)); 1/2, 1/sqrt(3); -1/2, 1/sqrt(3)], 'Space', 'vector');

% Re-define constitutive properties
v = 0;
E = 504;
D = E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]; % constitutive matrix

% Re-define Short-hand for shape function derivatives
B = @(xi1, xi2) elem.shape_deriv(xi1,xi2);

% Compute element stiffness matrix (compare to Eq. 24.37)
Ke = zeros(elem.n_dof);
for i = 1:length(qp);
    Ke = Ke + W(i)*B(qp(i,1),qp(i,2))'*D*B(qp(i,1),qp(i,2))*elem.detJ(qp(i,1),qp(i,2));
end
Ke


