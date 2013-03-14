function work


mesh = mFEM.Mesh();
mesh.grid('Quad4',0,1,0,1,10,5);
% mesh.addBoundary(1);
% mesh.update();
toc;
