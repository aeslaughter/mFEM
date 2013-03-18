function work

n = 13;
tic;
mesh = mFEM.Mesh();
mesh.grid('Quad4',0,1,0,1,n,n);
mesh.addBoundary(1);
mesh.update();
toc;
