function test_FunctionKernelRegistry
%TEST_FUNCTIONKERNELREGISTRY Tests the ConstantKernelRegistry class

% Call the Test class for this file
test = Test(mfilename('fullfile'));

% Run the code that is being tested
out = test.func(@code_to_test);

% Evaluate the results
test.compare(out(1), 0.25, 'String function added to registry');
test.compare(out(2), 70, 'String function with global constant registry');
test.compare(out(3), 120, 'String function with local constant registry');
test.compare(out(4), 520, 'String function with global constant registry, change value and re-evaluated');
test.compare(out(5), 120, 'Re-run #3, nothing changed in local constant registry');
end

function out = code_to_test

    constReg = mFEM.registry.ConstantKernelRegistry('-disableWarnings');
    constReg.add('a',5);
    constReg2 = mFEM.registry.ConstantKernelRegistry();
    constReg2.add('a',10);

    reg = mFEM.registry.FunctionKernelRegistry('Constants',constReg,'-disableWarnings');
    
    reg.add('func1','x(1)^2');
    reg.add('func2', 't*(x(1)+a)');
    reg.add('func3', 't*(x(1)+a)', 'constants', constReg2);

    elem = mFEM.elements.Line2(1,[0;1]);
    out(1) = reg.kernels(1).eval(elem,0.5);  % 0.5^2 = 0.25
    out(2) = reg.kernels(2).eval(elem,2,10); % 10*(2+5) = 70    
    out(3) = reg.kernels(3).eval(elem,2,10); % 10*(2+10) = 120
    
    constReg.add('a',50);
    out(4) = reg.kernels(2).eval(elem,2,10); % 10*(2+50) = 520
    out(5) = reg.kernels(3).eval(elem,2,10); % 10*(2+10) = 120

end