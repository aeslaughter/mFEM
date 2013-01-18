function test_Registry
%TEST_REGISTRY Tests the Registry class

% Call the Test class for this file
clear classes;
test = Test(mfilename('fullfile'));

% Run the code that is being tested
out = test.func(@code_to_test);

% Evaluate the results
test.compare(out{1},0.5, 'Function, with constant evaluated from registry');

end

function out = code_to_test
    reg = mFEM.registry.Registry();

    reg.add_constant('k',2)
    reg.add_function('fcn','k*x(1)^2');
    
    elem = mFEM.elements.Line2(1,[0;1]);
    out{1} = reg.functions.kernels(1).eval(elem,0.5);
end