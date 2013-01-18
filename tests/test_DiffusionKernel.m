function test_DiffusionKernel
%TEST_MATRIXKERNEL Tests the ConstantKernel class

% Call the Test class for this file
test = Test(mfilename('fullfile'));

% Run the code that is being tested
out = test.func(@code_to_test);

% Evaluate the results
M1 = [5.3125,-0.625,-4.6875; -0.625, 1.25, -0.625; -4.6875, -0.625, 5.3125];
M2 = [4.75961538461539 -3.50961538461538 -2.98076923076923 1.73076923076923;-3.50961538461538 4.13461538461539 1.73076923076923 -2.35576923076923;-2.98076923076923 1.73076923076923 6.53846153846154 -5.28846153846154;1.73076923076923 -2.35576923076923 -5.28846153846154 5.91346153846154]; 
test.compare(out{1}, M1, 'B''*D*B on Tri3, Fish, 2007, p. 194');
test.compare(out{2}, M2, 'B''*D*B on Quad4, Fish, 2007, p. 199', 'Tol', 10^-15);
end

function out = code_to_test

    % Define a single element for testing
    mesh = mFEM.FEmesh('Element','Tri3','time',false);
    mesh.add_element([0,0; 2,0.5; 0,1]);
    mesh.init();

    % TEST 1
    % K(1) from Fish (2007), p. 194 (Example 8.1)
    kern = mFEM.kernels.Diffusion(mesh,'D', 5);
    out{1} = kern.eval(mesh.element(1),[0.5,0.5]);   
    
    % TEST 2
    % K from Fish (2007), p. 199 (Example 8.2)
    mesh = mFEM.FEmesh('Element','Quad4','time',false);
    mesh.add_element([0,1; 0,0; 2,0.5; 2,1]);
    mesh.init();
    kern = mFEM.kernels.Diffusion(mesh,'D', 5);
    elem = mesh.element(1);
    [qp,W] = elem.quad.rules('-cell');
    Ke = zeros(elem.n_dof,elem.n_dof);
    for i = 1:length(qp)
        Ke = Ke + W(i)*kern.eval(elem,qp{i})*elem.detJ(qp{i});
    end
    out{2} = Ke;
end



