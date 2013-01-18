function test_ConstantKernelRegistry
%TEST_CONSTANTKERNELREGISTRY Tests the ConstantKernelRegistry class

% Call the Test class for this file
test = Test(mfilename('fullfile'));

% Run the code that is being tested
out = test.func(@code_to_test);

% Evaluate the results
test.compare(out.kernels(1).eval(),3, 'a = 3');
test.compare(out.kernels(2).eval(),2, 'b = 2');
test.compare(out.kernels(3).eval(),10, 'c = 10');

end

function varargout = code_to_test
    reg = mFEM.registry.base.ConstantKernelRegistry('-disableWarnings');

    reg.add('a','1');
    reg.add('b', 2);
    reg.add('a', 3);
    reg.add('c', '2*a');
    reg.add('c','4','-add');
    varargout{1} = reg;
end