function test

clear;
import mFEM.*

mesh = FEmesh('Element','Quad4');
mesh.grid(0,1,0,1,1,1);
mesh.init();
% mesh.plot()

sys = System(mesh);

sys.add_constant('a',0.1);

% sys.kernels{1}.init(sys.kernels);

sys.add_constant('b','2*a');

sys.kernels{1}.init();

sys.kernels{:}
return;

sys.add_constant('c','52');
sys.add_constant('d','a + b*c');
sys.add_constant('b',10);

sys.add_matrix('K','B''*a*B');
sys.add_matrix('K2', 'K + B''*b*B');
sys.add_matrix('K','B''*b*B');
%sys.add('Constant',...,'Const',...,'C',...,'matrix',...,'mat',...,'m',....,'vector,.......

sys.kernels{:}
