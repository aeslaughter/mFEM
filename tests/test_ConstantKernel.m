function test_ConstantKernel
%TEST_CONSTANTKERNEL Tests the ConstantKernel class

% Call the Test class for this file
test = Test(mfilename('fullfile'));

% Run the code that is being tested
out = test.func(@code_to_test);

% Evaluate the results
test.compare(out{1}.eval(),0.1, 'k = 0.1');
test.compare(out{2}.eval(),0.2, 'D = 0.2');

end

function c = code_to_test
    c{1} = mFEM.kernels.base.ConstantKernel('k',0.1);
    c{2} = mFEM.kernels.base.ConstantKernel('D','0.2');
end



