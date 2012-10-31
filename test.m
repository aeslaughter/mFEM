function test
import mFEM.*

mesh = FEmesh();
mesh.grid('Tri6',0,1,0,1,1,1);
mesh.init();


% 
%  mesh.element(1).quad.rules()
T = 1:mesh.n_dof;
mesh.plot(T,'ShowNodes',true);