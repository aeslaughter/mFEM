function T = test_FunctionKernelRegistry
%TEST_FUNCTIONKERNELREGISTRY Tests the ConstantKernelRegistry class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Create an element
elem = mFEM.elements.Line2(1,[0;1]);

% Create two constant registry objects
constReg = mFEM.registry.ConstantRegistry('-disableWarnings');
constReg.add('a',5);
constReg2 = mFEM.registry.ConstantRegistry();
constReg2.add('a',10);

% Create a function registry object
reg = mFEM.registry.FuncRegistry('ConstantRegistry',constReg,'-disableWarnings');

% Perform tests
reg.add('func1', 'x(1)^2', 'func2', 't*(x(1)+a)');
val = reg.get('func1', elem, 0.5);  % 0.5^2 = 0.25
T.compare(val, 0.25, 'String function added to registry');

%reg.add('func2', 't*(x(1)+a)');
val = reg.get('func2',elem,2,10); % 10*(2+5) = 70    
T.compare(val, 70, 'String function using global constant registry');

reg.add('func3', 't*(x(1)+a)','ConstantRegistry',constReg2);
val = reg.get('func3',elem,2,10); % 10*(2+10) = 120
T.compare(val, 120, 'String function with local constant registry');

constReg.add('a',50);
val = reg.get('func2',elem,2,10); % 10*(2+50) = 520
T.compare(val, 520, 'String function with global constant registry, change value and re-evaluated');

% Evaluate the results
val = reg.get('func3',elem,2,10); % 10*(2+10) = 120
T.compare(val, 120, 'Validate that nothing changed in local constant registry');




