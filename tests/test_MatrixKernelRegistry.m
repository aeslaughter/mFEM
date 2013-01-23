function T = test_MatrixKernelRegistry
%TEST_MATRIXKERNEL Tests the ConstantKernel class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Create a mesh (single Tri3 element);
mesh = mFEM.FEmesh('Element','Tri3','time',false);
mesh.add_element([0,0; 2,0.5; 0,1]);
mesh.init();

% Create two Diffusion kernels (each 1/2 of what a single should be) for
% computing K(1) from Fish (2007), p. 194 (Example 8.1)
kern1 = mFEM.kernels.Diffusion(mesh,'D', 2.5);
kern2 = mFEM.kernels.Diffusion(mesh,'D', 2.5);

% Create a Kernel Registry and add the two kernels
reg = mFEM.registry.MatrixKernelRegistry(mesh);
reg.add('K', kern1);
reg.add('K', kern2);

% Assemble the 'K' matrix and test that it is correct
Kexact = [5.3125,-0.625,-4.6875; -0.625, 1.25, -0.625; -4.6875, -0.625, 5.3125];
Kcalc = reg.assemble('K');
T.compare(Kexact, Kcalc, 'B''*D*B on Tri3, Fish, 2007, p. 194');

% Test the handle and associated matrix are identical for the two kernels
T.compare(kern1.get(),kern2.get(), 'Matrix numbers are identical');
T.compare(kern1.value, kern2.value, 'Matrix handles are identical');
