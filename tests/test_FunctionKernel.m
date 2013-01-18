function test_FunctionKernel
%TEST_FUNCTIONKERNEL Tests the ConstantKernel class

% Call the Test class for this file
test = Test(mfilename('fullfile'));

% Run the code that is being tested
out = test.func(@code_to_test);

% Evaluate the results
test.compare(out{1}, 8, 'Inline function');
test.compare(out{2}, 24, 'Text function');
test.compare(out{3}, 27, 'Externel, element function');
end

function out = code_to_test
    elem = mFEM.elements.Line2(1,[0;3]);

    kern = mFEM.kernels.base.FunctionKernel('func', @(elem,x,t) 2*x^2);
    out{1} = kern.eval(elem,2,0);
    
    kern = mFEM.kernels.base.FunctionKernel('func', '3*x^3');
    out{2} = kern.eval(elem,2,0);
    
    fcn = @(elem,x,t) elem.size()*t*x;
    kern = mFEM.kernels.base.FunctionKernel('func', fcn);
    out{3} = kern.eval(elem,3,3);
end



