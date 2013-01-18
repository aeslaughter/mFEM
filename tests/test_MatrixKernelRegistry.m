function test_MatrixKernelRegistry
%TEST_MATRIXKERNEL Tests the ConstantKernel class

% Call the Test class for this file
test = Test(mfilename('fullfile'));

% Run the code that is being tested
out = test.func(@code_to_test);

% Evaluate the results
% test.compare(out{1}, 8, 'Inline function');
% test.compare(out{2}, 24, 'Text function');
% test.compare(out{3}, 27, 'Externel, element function');
end

function out = code_to_test
    % Define a single element for testing
    mesh = mFEM.FEmesh('Element','Tri3','time',false);
    mesh.add_element([0,0; 2,0.5; 0,1]);
    mesh.init();

    % TEST 1
    % K(1) from Fish (2007), p. 194 (Example 8.1)
    kern1 = mFEM.kernels.Diffusion(mesh,'D', 5);
    kern2 = mFEM.kernels.Diffusion(mesh,'D', 5);
    kern3 = mFEM.kernels.Diffusion(mesh,'D', 5);
    
    reg = mFEM.registry.MatrixKernelRegistry(mesh);
    reg.add('K', kern1);
    reg.add('K2', kern2);
    reg.add('K', kern3);
    
    reg.locate('K')
    
    
    out(1) = 0;
%     out{1} = kern.eval(elem,2,0);
%     
%     kern = mFEM.kernels.base.FunctionKernel('func', '3*x^3');
%     out{2} = kern.eval(elem,2,0);
%     
%     fcn = @(elem,x,t) elem.size()*t*x;
%     kern = mFEM.kernels.base.FunctionKernel('func', fcn);
%     out{3} = kern.eval(elem,3,3);
end



