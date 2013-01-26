function T = test_example1

% Create Test class
T = mFEM.Test();

% Compute exact solution
x = 0:0.2:4;
yex = -12.5*x.^2 + 97.5*x;

% Solve example1a
y = example1a(length(x)-1);
T.compare(y,yex','Example1a: Temperature solution correct', 'Tol', 10^-12);

% Solve example1b
y = example1b('N',length(x)-1);
T.compare(y,yex','Example1b: Temperature solution correct', 'Tol', 10^-12);

% Solve example1c
y = example1c(length(x)-1);
T.compare(y,yex','Example1c: Temperature solution correct', 'Tol', 10^-12);
h = findobj('Type','Figure','-and','IntegerHandle','on');
close(h);