function T = test_example5

% Initilize testing class
T = mFEM.Test();


% Run example 5a
[Ma,Ka,fa,Ta] = example5a('-debug');
performTests(T, 'Example5a', Ma, Ka, fa, Ta);

[Mb,Kb,fb,Tb] = example5b('-debug');
performTests(T, 'Example5b', Mb, Kb, fb, Tb);

function performTests(x, name, M, K, f, T)

% "Exact solution," taken from known working set, see also Bhatti 2005.
Tex = [50 145.031526140504 160.751903085136 175.901546899857 180.149547931426 182.806447442623 183.754905085443 184.249105278178 184.447125625046 184.542048394261 184.58224523089;50 78.2358744154568 123.605675519846 133.756244557911 141.356159096001 143.796177047445 145.171916188397 145.695831329464 145.956281108166 146.063924266308 146.11440560176;300 300 300 300 300 300 300 300 300 300 300;300 300 300 300 300 300 300 300 300 300 300];
Mex = sparse([1;2;3;4;1;2;3;1;2;3;4;1;3;4], [1;1;1;1;2;2;2;3;3;3;3;4;4;4], [128;42.6666666666667;64;21.3333333333333;42.6666666666667;85.3333333333333;42.6666666666667;64;42.6666666666667;128;21.3333333333333;21.3333333333333;21.3333333333333;42.6666666666667], 4, 4);
Kex = sparse([1;2;3;4;1;2;3;1;2;3;4;1;3;4], [1;1;1;1;2;2;2;3;3;3;3;4;4;4], [7;-2;1.5;-4.5;-2;4.75;-0.75;1.5;-0.75;2.25;-3;-4.5;-3;7.5], 4, 4);
fex = [100;100;0;0];

tol = 10^-13;
x.compare(Tex(1,:), T(1,:), [name,': First time-step temperature.'], 'Tol', 10^-12);
x.compare(Tex(end,:), T(end,:), [name, ': Last time-step temperature.'], 'Tol', tol);
x.compare(Mex, M, [name,': Mass matrix'], 'Tol', tol);
x.compare(Kex, K, [name,': Stiffness matrix'], 'Tol', tol);
x.compare(fex, f, [name,': Force vector']);