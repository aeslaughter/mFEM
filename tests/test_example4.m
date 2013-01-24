function test_example4

T = mFEM.Test();

% "Exact" stress and strain; values from a known working solution
strain_exact = [8.82024842022898e-07 -3.61335341871653e-07 6.6511117087177e-07 -1.17086818346395e-06;-6.27568703391555e-08 -6.27568703391555e-08 -3.45843535505952e-07 -3.45843535505952e-07;-4.0260763424007e-06 -3.94034886325541e-06 9.46655214265826e-08 2.21253041156388e-07];
stress_exact = [28.4570697006973 -12.5328264716505 18.5063113259335 -42.020480575244;6.65441480003452 -5.64255405166983 -4.82341266739851 -22.9814502377518;-46.4547270277004 -45.4655638067932 1.09229447799903 2.55291970565063];

[stress, strain] = example4a('-display');
tol = 10*eps;
T.compare(stress, stress_exact, 'Example4a: Stress', 'Tol', tol);
T.compare(strain, strain_exact, 'Example4a: Strain', 'Tol', tol);

[stress, strain] = example4b('-display');
T.compare(stress, stress_exact, 'Example4b: Stress', 'Tol', tol);
T.compare(strain, strain_exact, 'Example4b: Strain', 'Tol', tol);
