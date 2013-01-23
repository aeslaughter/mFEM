function T = test_Diffusion
%TEST_Diffusion Tests the Diffusion kernel

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Define a single element mesh (Tri3)
mesh = mFEM.FEmesh('Element','Tri3','time',false);
mesh.add_element([0,0; 2,0.5; 0,1]);
mesh.init();

% K(1) from Fish (2007), p. 194 (Example 8.1)
kern = mFEM.kernels.Diffusion(mesh,'D', 5);
Kexact = [5.3125,-0.625,-4.6875; -0.625, 1.25, -0.625; -4.6875, -0.625, 5.3125];
Kcalc = kern.eval(mesh.element(1),[0.5,0.5]);   
T.compare(Kexact, Kcalc, 'B''*D*B on Tri3, Fish, 2007, p. 194');

% Define a single element mesh (Quad4)
mesh = mFEM.FEmesh('Element','Quad4','time',false);
mesh.add_element([0,1; 0,0; 2,0.5; 2,1]);
mesh.init();

% K from Fish (2007), p. 199 (Example 8.2)
kern = mFEM.kernels.Diffusion(mesh,'D', 5);
Kexact = [4.75961538461539 -3.50961538461538 -2.98076923076923 1.73076923076923;-3.50961538461538 4.13461538461539 1.73076923076923 -2.35576923076923;-2.98076923076923 1.73076923076923 6.53846153846154 -5.28846153846154;1.73076923076923 -2.35576923076923 -5.28846153846154 5.91346153846154]; 
Kcalc = kern.assemble('-zero');
T.compare(Kexact, Kcalc, 'B''*D*B on Quad4, Fish, 2007, p. 199', 'Tol', 10^-15);

% Test that the matrix was cleared
Kcalc = kern.assemble();
T.compare(Kexact, Kcalc, 'Matrix properly cleared', 'Tol', 10^-15);

