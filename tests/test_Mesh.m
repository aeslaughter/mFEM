function test_Mesh

mesh = mFEM.Mesh();
mesh.grid('Quad8',0,1,0,1,2,1);