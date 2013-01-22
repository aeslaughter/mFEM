function T = test_Func
%TEST_FUNC Tests the Func class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Create an element
elem = mFEM.elements.Line2(1,[0;3]);

kern = mFEM.kernels.Func('func', @(elem,x,t) 2*x^2);
T.compare(kern.eval(elem,2,0), 8, 'Direct function handle input');

kern = mFEM.kernels.Func('func', '3*x^3');
T.compare(kern.eval(elem,2,0), 24, 'Function string input');

fcn = @(elem,x,t) elem.size()*t*x;
kern = mFEM.kernels.Func('func', fcn);
T.compare(kern.eval(elem,3,3), 27, 'Externel function handle input');

% constReg = mFEM.registry.ConstantKernelRegistry();
% constReg.add('a',5);
% kern = mFEM.kernels.base.Func('func', 't*(x(1)+a)', 'constants', constReg);
% T.compare(kern.eval(elem,2,10), 70, 'String function with constant from registry');
