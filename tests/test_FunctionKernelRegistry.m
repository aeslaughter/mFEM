function test_FunctionKernelRegistry
%TEST_FUNCTIONKERNELREGISTRY Tests the ConstantKernelRegistry class

% Call the Test class for this file
test = Test(mfilename('fullfile'));

% Run the code that is being tested
out = test.func(@code_to_test);

% Evaluate the results
test.compare(out{1},0.25, 'Function added to registry');

end

function out = code_to_test
    reg = mFEM.registry.base.KernelRegistry('Type','function','-disableWarnings');

    reg.add('f1','x(1)^2');
    
    elem = mFEM.elements.Line2(1,[0;1]);
    out{1} = reg.kernels(1).eval(elem,0.5);
end