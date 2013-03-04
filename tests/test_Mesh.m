function test_Mesh

tic;
mesh = mFEM.Mesh();
% profile on;
% mFEM.elements.Quad4.buildMaps(0,1,0,1,5000,1000);
% profile viewer;
mesh.grid('Quad4',0,1,0,1,10,10);
toc