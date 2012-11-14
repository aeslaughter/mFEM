function test

clear;
import mFEM.*

mesh = FEmesh();
mesh.add_element('Truss',[0,0; 8.66,0]);
mesh.init();

elem = mesh.element(1);
N = elem.shape();
T = elem.transformation()

L = elem.size()
E = 1e6;
A = 0.01;

A*E/L*N'*N




return;



B = @(xi) elem.shape(xi);



[qp,W] = elem.quad.rules(); 
K = zeros(elem.n_dof*2)
for i = 1:length(qp);
    K = K + W(i)*A*E/L*B(qp(i))'*B(qp(i))*elem.detJ(qp(i));
end
K
% sys = System(mesh);
% sys.add_matrix('K','B''*B');
% 
% K = sys.assemble('K'); full(K)













% mesh = FEmesh();
% mesh.grid('Quad4',0,1,0,1,1,1);
% mesh.init();

% sys = System(mesh);
% sys.add_constant('k', 10);
% sys.add_matrix('K','B''*k*B');
% %sys.add_matrix('K','B''*k*B');

% K = sys.assemble('K'); full(K)

% func.fhandle{1} = @(x,t) test_fcn(x,t);
% 
% Fstr = '@(elem, x, xi, eta) elem.shape(xi,eta)''*feval(func.fhandle{1},x,1)*elem.shape(xi,eta)'
% F = str2func(Fstr);
% 
% M = Matrix(mesh);
% for e = 1:mesh.n_elements;
%     elem = mesh.element(e);
%     [qp,w] = elem.quad.rules();
%     
%     Me = zeros(elem.n_dof);
%     for i = 1:length(qp);
%         for j = 1:length(qp);
%             x = elem.get_position(qp(i),qp(j));
%             Me = Me + w(i)*w(j)*F(elem,x,qp(i),qp(j))*elem.detJ(qp(i),qp(j));
%         end
%     end
%     
%     M.add_matrix(Me, elem.get_dof());
%     
% end
% 
% full(M.init())

% sys = System(mesh);
% sys.add_constant('c', 1);
% sys.add_function('k', @(x,t) test_fcn(x,t));
% sys.time = 2;
% sys.add_matrix('M', 'N''*c*k*N');
% 
% M = sys.assemble('M');
% % full(M)
% 
% function output = test_fcn(x,t)
%     output = 2;
%     %output = 2*x(1)*x(2)*t;

