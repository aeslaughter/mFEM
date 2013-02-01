function T = test_System
%TEST_SYSTEM Tests the ConstantKernel class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Define a single element mesh (Tri3)
mesh = mFEM.FEmesh('time', false);
mesh.addElement('Tri3',[0,0; 2,0.5; 0,1]);
mesh.init();

% K(1) from Fish (2007), p. 194 (Example 8.1)
sys = mFEM.System(mesh);

% Add Constants
sys.addConstant('D', 15);
T.compare(sys.get('D'), 15, 'Create and retrieve a constant.');

% Add to a consant
sys.addConstant('D', 5,'-add');
T.compare(sys.get('D'), 20, 'Add to a a constant.');

% Overwrite constant
sys.addConstant('D', 5);
T.compare(sys.get('D'), 5, 'Overwrite a constant.');

% Add vector



%
sys.addMatrix('K', 'B''*D*B');
Kcalc = sys.assemble('K');
Kexact = [5.3125,-0.625,-4.6875; -0.625, 1.25, -0.625; -4.6875, -0.625, 5.3125];
T.compare(Kexact, Kcalc, 'B''*D*B on Tri3, Fish, 2007, p. 194');


