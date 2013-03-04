function T = test_example5

% Initilize testing class
T = mFEM.Test();


% Run example 5a
[Ma,Ka,fa,Ta] = example5a('-debug');
performTests(T, 'Example5a', Ma, Ka, fa, Ta);

[Mb,Kb,fb,Tb] = example5b('-debug');
performTests(T, 'Example5b', Mb, Kb, fb, Tb);

function performTests(x, name, M, K, f, T)

% "Exact solution" taken from known working set, Mex and Kex taken from  Bhatti 2005.
Tex = [50 145.031526140504 160.751903085136 175.901546899857 180.149547931426 182.806447442623 183.754905085443 184.249105278178 184.447125625046 184.542048394261 184.58224523089;50 78.2358744154568 123.605675519846 133.756244557911 141.356159096001 143.796177047445 145.171916188397 145.695831329464 145.956281108166 146.063924266308 146.11440560176;300 300 300 300 300 300 300 300 300 300 300;300 300 300 300 300 300 300 300 300 300 300];
Mex = [128,64/3,128/3,64; 64/3,128/3,0,64/3; 128/3,0,256/3,128/3; 64,64/3,128/3,128];
Kex = [22/3,-9/2,-7/3,3/2; -9/2,15/2,0,-3; -7/3,0,61/12,-3/4; 3/2,-3,-3/4,9/4];
fex = [100;0;100;0];
ix = [1,4,2,3];

tol = 10^-12;
x.compare(Tex(1,:), T(1,:), [name,': First time-step temperature.'], 'Tol', tol);
x.compare(Tex(end,:), T(end,:), [name, ': Last time-step temperature.'], 'Tol', tol);
x.compare(Mex, M(ix,ix), [name,': Mass matrix'], 'Tol', tol);
x.compare(Kex, K(ix,ix), [name,': Stiffness matrix'], 'Tol', tol);
x.compare(fex, f(ix,:), [name,': Force vector']);