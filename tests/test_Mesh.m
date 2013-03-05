function test_Mesh

mesh = mFEM.Mesh();
mesh.grid('Quad4',0,1,0,1,20,20);
mesh.plot([]);



delete(mesh);

% 344x300 fails: pkg = 80392 bytes (343x300 works)
