function T = test_example3

% Initilize testing class
T = mFEM.Test();

% Run example 5a
[Ta,qa] = example3a('-debug');
performTests(T, 'Example3a', Ta, qa);

[Tb,qb] = example3b('-debug');
performTests(T, 'Example3b', Tb, qb);

function performTests(x, name, T, q)

% "Exact solution," taken from known working set, see also Bhatti 2005.
Tex = [0;0;0;-3.04390243902439];
qex = [5.81163115918566 0.899062469187655 4.95459606739771 -2.29943603723445;3.59624987675062 3.59624987675062 19.8183842695908 19.8183842695908];

tol = 10^-13;
x.compare(Tex, T, [name,': temperature.'], 'Tol', 10^-12);
x.compare(qex, q, [name, ': flux.'], 'Tol', tol);