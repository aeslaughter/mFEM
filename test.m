function test

clear;
import mFEM.*

mesh = FEmesh();
mesh.grid('Quad4',0,1,0,1,1,1);
mesh.init();

sys = System(mesh);
sys.add_constant('k', 10);
sys.add_matrix('K','B''*k*B');
%sys.add_matrix('K','B''*k*B');

K = sys.assemble('K'); full(K)

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

function output = test_fcn(x,t)
    output = 2;
    %output = 2*x(1)*x(2)*t;

