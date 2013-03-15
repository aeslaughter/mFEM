function T = test_example13

T = mFEM.Test();
[~,r] = example13a('-debug');
rex = [-16;23;33]*1000;
T.compare(r,rex, 'Example13a: Reactions', 'Tol', 10^-8);

[~,r] = example13b('-debug');
T.compare(r,rex, 'Example13b: Reactions', 'Tol', 10^-8);
