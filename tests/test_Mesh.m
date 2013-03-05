function test_Mesh

mesh = mFEM.Mesh();
mesh.grid('Quad4',0,1,0,1,1000,1000);
