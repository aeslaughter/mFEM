function T = test_example2

% Initilize testing class
T = mFEM.Test();

% Run example 5a
[Ta,qa] = example2a('-debug');
performTests(T, 'Example2a', Ta, qa);

[Tb,qb] = example2b('-debug');
performTests(T, 'Example2b', Tb, qb);

function performTests(x, name, T, q)

% "Exact solution," taken from known working set, see also Bhatti 2005.
Tex = [0;0;0;-1.78823529411765];
qex = [0 4.47058823529412;0 17.8823529411765];

tol = 10^-13;
x.compare(Tex, T, [name,': temperature.'], 'Tol', 10^-12);
x.compare(qex, q, [name, ': flux.'], 'Tol', tol);