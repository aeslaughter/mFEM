function T = test_example6

% Initilize testing class
T = mFEM.Test();

% Run example 6a
[x,y,t,Temp] = example6a('-debug', 'N', 8);
performTests(T, 'Example6a', x, y, t, Temp);

[x,y,t,Temp]  = example6b('-debug', 'N', 8);
performTests(T, 'Example6b', x, y, t, Temp);

function performTests(C, name, x, y, t, T)

% Exact solution
Tex = @(x,y,t) exp(-t)*sin(pi*x).*sin(pi*y);
C.compare(Tex(x,y,t), T, [name,': Temperature solution.'], 'Tol', 10^-13);
