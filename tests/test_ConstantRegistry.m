function T = test_ConstantRegistry
%TEST_CONSTANTREGISTRY Tests the ConstantKernelRegistry class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Create a registry
reg = mFEM.registry.ConstantRegistry('-disableWarnings');

reg.add('a','2');
T.compare(reg.get('a'), 2, 'Adding kernel, text input');

reg.add('b', 5);
T.compare(reg.get('b'), 5, 'Adding kernel, numeric input');

reg.add('c', '2*a'); 
T.compare(reg.get('c'), 4, 'Adding kernel, function of another constant');

reg.add('b', 10);
T.compare(reg.get('b'), 10, 'Adding kernel, replace existint value');

reg.add('d', 32, 'e', 43');
T.compare(reg.get('d'), 32, 'Adding kernel, first input of multiple inputs');
T.compare(reg.get('e'), 43, 'Adding kernel, second input of multiple inputs');

reg.add('q', 40, '-add');
T.compare(reg.get('q'), 40, 'Add flag, not repated value');

reg.add('d',10, '-add'); % 32 + 10 = 42
T.compare(reg.get('d'),42, 'Added numeric value to existing constant');

reg.add('d','2*b + c','-add'); % 42 + 2*10 + 2*2 = 66
T.compare(reg.get('d'),66, 'Added string function to existing constant');
end