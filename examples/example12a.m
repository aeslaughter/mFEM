%EXAMPLE12a A truss example

function example12a


import mFEM.*

mesh = FEmesh('Space','Truss');
mesh.add_element('Truss',[0,0; 8.66,0]);
mesh.add_element('Truss',[8.66,0; 0,5]);
mesh.add_element('Truss',[0,5; 0,0]);
mesh.init();

E = 1e6;
A = 0.01;
K = Matrix(mesh);
for e = 1:mesh.n_elements;
    elem = mesh.element(e)
    
    L = elem.size();
    T = elem.transformation();
    N = elem.shape()
    Ke = A*E/L*T'*(N'*N)*T
    elem.get_dof()
    K.add_matrix(Ke,elem.get_dof());
    
    
end

K = K.init(); full(K)

idx = [3,4,5,6,1,2];

full(K(idx,idx))





