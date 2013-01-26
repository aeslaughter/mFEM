function T = test_example7

T = mFEM.Test();
[u,uex] = example7b('-debug');
T.compare(u,uex, 'Displacement', 'tol', 10^-2);