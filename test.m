function test

clear;
import mFEM.*

mesh = FEmesh();
mesh.grid('Quad4',0,1,0,1,3,3);
mesh.init();

M = Matrix(mesh);

fhandle = @(x,y,t) fcn(x,y,t);

for e = 1:mesh.n_elements;
   elem = mesh.element(e);
   [qp,w] = elem.quad.rules();
   Me = zeros(elem.n_dof);
   for i = 1:length(qp);
       for j = 1:length(qp);
           [x,y] = elem.get_position(qp(i),qp(j));
           Me = Me + elem.shape(qp(i),qp(j))'*elem.shape(qp(i),qp(j))*feval(fhandle,x,y,1)*elem.detJ(qp(i),qp(j));
       end
   end
   M.add_matrix(Me,elem.get_dof());
end

full(M.init())

function output = fcn(x,y,t)
    output = 2*x*y*t;

