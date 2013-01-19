function T = test_ConstantKernelRegistry
%TEST_CONSTANTKERNELREGISTRY Tests the ConstantKernelRegistry class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Create a registry
reg = mFEM.registry.ConstantKernelRegistry('-disableWarnings');

reg.add('a','2');
T.compare(reg.kernels(1).eval(), 2, 'Adding kernel, text input');

reg.add('b', 5);
T.compare(reg.kernels(2).eval(), 2, 'Adding kernel, numeric input');

reg.add('c', '2*a'); 
T.compare(reg.kernels(3).eval(), 2, 'Adding kernel, function of another constant');

reg.add('b', 10);
T.compare(reg.kernels(2).eval(), 10, 'Adding kernel, replace existint value');

reg.add('d', 32, 'e', 43');
T.compare(reg.kernels(4).eval(), 10, 'Adding kernel, first input of multiple inputs');
T.compare(reg.kernels(5).eval(), 10, 'Adding kernel, second input of multiple inputs');

end