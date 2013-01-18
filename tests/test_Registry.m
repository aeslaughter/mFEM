function test_Registry
%TEST_REGISTRY Tests the Registry class

% Call the Test class for this file
clear all;
test = Test(mfilename('fullfile'));

% Run the code that is being tested
out = test.func(@code_to_test);

% Evaluate the results
test.compare(out{1},0.5, 'Function, with constant evaluated from registry');
test.compare(out{2},16, 'Function handle evaluated from registry');
test.compare(out{3},4, 'Constant evaluated from registry');

end

function out = code_to_test
    reg = mFEM.registry.Registry();

    reg.add_constant('k',2)
    reg.add_constant('D','2*k');
    reg.add_function('fcn','k*x(1)^2');
    reg.add_function('fcn2',@(elem,x,t) x(1)^4);
    
    elem = mFEM.elements.Line2(1,[0;1]);
    out{1} = reg.functions.kernels(1).eval(elem,0.5);
    out{2} = reg.functions.kernels(2).eval(elem,2);
    out{3} = reg.constants.kernels(2).eval();
    
    delete('reg');
end