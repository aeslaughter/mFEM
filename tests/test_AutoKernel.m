function T = test_AutoKernel(varargin)
%TEST_AUTOKERNEL Tests the ConstantKernel class

% Call the Test class for this file
T = mFEM.Test('Name','AutoKernel',varargin{:});

% Define a single element mesh (Tri3)
mesh = mFEM.Mesh();
mesh.createNode([0,0]);
mesh.createNode([2,0.5]);
mesh.createNode([0,1]);
mesh.createElement('Tri3',[1,2,3]);
mesh.init();

% K(1) from Fish (2007), p. 194 (Example 8.1)
kern = mFEM.kernels.AutoKernel(mesh,'K','B''*D*B', 'D', 5);
Kexact = [5.3125,-0.625,-4.6875; -0.625, 1.25, -0.625; -4.6875, -0.625, 5.3125];
Kcalc = kern.eval(mesh.element(1),[0.5,0.5],[]);   
T.compare(Kexact, Kcalc, 'B''*D*B on Tri3, Fish, 2007, p. 194');
return;
% Define a single element mesh (Quad4)
mesh2 = mFEM.FEmesh('time', false);
mesh2.addElement('Quad4',[0,1; 0,0; 2,0.5; 2,1]);
mesh2.init();

% K from Fish (2007), p. 199 (Example 8.2)
kern = mFEM.kernels.AutoKernel(mesh2,'K','B''*D*B', 'D', 5);
Kexact2 = [4.75961538461539 -3.50961538461538 -2.98076923076923 1.73076923076923;-3.50961538461538 4.13461538461539 1.73076923076923 -2.35576923076923;-2.98076923076923 1.73076923076923 6.53846153846154 -5.28846153846154;1.73076923076923 -2.35576923076923 -5.28846153846154 5.91346153846154]; 
Kcalc = kern.assemble('-zero');
%T.compare(Kexact2, Kcalc, 'B''*D*B on Quad4, Fish, 2007, p. 199', 'Tol', 10^-15);

% Test that the matrix was cleared
Kcalc = kern.assemble();
T.compare(Kexact2, Kcalc, 'Matrix properly cleared', 'Tol', 10^-15);

% Test that constants registry is added properly
const = mFEM.registry.ConstantRegistry();
const.add('D', 5);
kern = mFEM.kernels.AutoKernel(mesh,'K', 'B''*D*B', 'ConstantRegistry', const);
Kcalc = kern.eval(mesh.element(1), [0.5,0.5],[]);
T.compare(Kexact, Kcalc, 'ConstantRegistry linking correctly');

% Test direct assemble
mesh = mFEM.FEmesh();
mesh.addElement('Truss',[0,0; 1,1]);
mesh.init();
kern = mFEM.kernels.AutoKernel(mesh,'K','Ke');
Kexact = 1/2*[1,1,-1,-1; 1,1,-1,-1; -1,-1,1,1; -1,-1,1,1];
Kcalc = kern.eval(mesh.element(1), [], []);
T.compare(Kexact, Kcalc, 'Direct stiffness assembly', 'Tol', 10^-14);


