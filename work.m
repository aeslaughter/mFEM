function work

tic;
mesh = mFEM.Mesh();
mesh.grid('Quad4',0,1,0,1,1000,200);
% mesh.addBoundary(1);
% mesh.update();
toc;
