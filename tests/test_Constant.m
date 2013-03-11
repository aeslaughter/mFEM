function T = test_Constant(varargin)
%TEST_CONSTANT Tests the Constant class

% Call the Test class for this file
T = mFEM.Test('Name','Constant',varargin{:});

% Run the code that is being tested
kern1 = mFEM.kernels.Constant('k', 0.1);
T.compare(kern1.eval(),0.1, 'Constant creation, numeric input');

% Evaluate the results
kern2 = mFEM.kernels.Constant('D','0.2');
T.compare(kern2.eval(),0.2, 'Constant creation, text input');

% Done
delete(T);
end




