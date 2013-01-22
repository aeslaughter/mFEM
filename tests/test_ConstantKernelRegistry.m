function T = test_ConstantKernelRegistry
%TEST_CONSTANTKERNELREGISTRY Tests the ConstantKernelRegistry class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Create a registry
reg = mFEM.registry.ConstantKernelRegistry('-disableWarnings');

reg.addConstant('a','2');
T.compare(reg.const(1).eval(), 2, 'Adding kernel, text input');

reg.addConstant('b', 5);
T.compare(reg.const(2).eval(), 5, 'Adding kernel, numeric input');

reg.addConstant('c', '2*a'); 
T.compare(reg.const(3).eval(), 4, 'Adding kernel, function of another constant');

reg.addConstant('b', 10);
T.compare(reg.const(2).eval(), 10, 'Adding kernel, replace existint value');

reg.addConstant('d', 32, 'e', 43');
T.compare(reg.const(4).eval(), 32, 'Adding kernel, first input of multiple inputs');
T.compare(reg.const(5).eval(), 43, 'Adding kernel, second input of multiple inputs');

end