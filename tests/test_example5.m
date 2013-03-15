function T = test_example5(varargin)

% Initilize testing class
T = mFEM.Test('Name','Example 5',varargin{:});

% Run example 5a
[Ma,Ka,fa,Ta] = example5a('-debug');
performTests(T, 'Example5a', Ma, Ka, fa, Ta);

[Mb,Kb,fb,Tb] = example5b('-debug');
performTests(T, 'Example5b', Mb, Kb, fb, Tb);

function performTests(x, name, M, K, f, T)

% "Exact solution," taken from known working set, see also Bhatti 2005.
Tex = [50 143.559425028756 159.066944434047 174.932967660026 178.906058331933 181.732370840661 182.62515163195 183.147314074371 183.336518536752 183.435439600772 183.474351107119;300 300 300 300 300 300 300 300 300 300 300;50 80.3726317160157 126.314513716304 135.598164612113 143.557693881035 145.780741042317 147.221844480129 147.707204473682 147.976494223602 148.077881744897 148.129282370427;300 300 300 300 300 300 300 300 300 300 300];
Mex = [128 21.3333333333333 42.6666666666667 64;21.3333333333333 42.6666666666667 0 21.3333333333333;42.6666666666667 0 85.3333333333333 42.6666666666667;64 21.3333333333333 42.6666666666667 128];
Kex = [7.33333333333333 -4.5 -2.33333333333333 1.5;-4.5 7.5 0 -3;-2.33333333333333 0 5.08333333333333 -0.75;1.5 -3 -0.75 2.25];
fex = [100;0;100;0];

tol = 10^-13;
x.compare(Tex(1,:), T(1,:), [name,': First time-step temperature.'], 'Tol', 10^-12);
x.compare(Tex(end,:), T(end,:), [name, ': Last time-step temperature.'], 'Tol', tol);
x.compare(Mex, M, [name,': Mass matrix'], 'Tol', tol);
x.compare(Kex, K, [name,': Stiffness matrix'], 'Tol', tol);
x.compare(fex, f, [name,': Force vector']);