function T = test_ConstantKernel
%TEST_CONSTANTKERNEL Tests the ConstantKernel class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Run the code that is being tested
kern1 = mFEM.kernels.base.ConstantKernel('k', 0.1);
T.compare(kern1.eval(),0.1, 'Constant creation, numeric input');

% Evaluate the results
kern2 = mFEM.kernels.base.ConstantKernel('D','0.2');
T.compare(kern2.eval(),0.2, 'Constant creation, text input');

end




