function test_ConstantKernelRegistry
%TEST_CONSTANTKERNELREGISTRY Tests the ConstantKernelRegistry class

% Call the Test class for this file
test = Test(mfilename('fullfile'));

% Run the code that is being tested
out = test.func(@code_to_test);

% Evaluate the results
test.compare(out(1), 2, 'Text input');
test.compare(out(2), 5, 'Numeric input');
test.compare(out(3), 4, 'Function of another input');
test.compare(out(4), 7, 'Replace existing');
end

function out = code_to_test
    reg = mFEM.registry.base.KernelRegistry('Type','constant','-disableWarnings');
    %reg = mFEM.registry.base.ConstantKernelRegistry();

    reg.add('a','2');
    reg.add('b', 5);
    reg.add('c', '2*a'); 
    reg.add('d', 5);
    reg.add('d', 7);
    
    for i = 1:length(reg.kernels);
        out(i) = reg.kernels(i).eval();
    end
end