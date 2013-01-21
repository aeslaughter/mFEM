function T = test_System
%TEST_SYSTEM Tests the ConstantKernel class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Define a single element mesh (Tri3)
mesh = mFEM.FEmesh('Element','Tri3','time',false);
mesh.add_element([0,0; 2,0.5; 0,1]);
mesh.init();

% K(1) from Fish (2007), p. 194 (Example 8.1)
sys = mFEM.System(mesh);
sys.add_constant('D',5);
sys.add_matrix('K','B''*D*B');
Kcalc = sys.assemble('K');
Kexact = [5.3125,-0.625,-4.6875; -0.625, 1.25, -0.625; -4.6875, -0.625, 5.3125];
T.compare(Kexact, Kcalc, 'B''*D*B on Tri3, Fish, 2007, p. 194');