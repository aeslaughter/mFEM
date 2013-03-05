function test_Mesh

mesh = mFEM.Mesh('-time');
mesh.grid('Quad4',0,2,0,2,2,2);
% mesh.plot([]);


delete(mesh);

% 344x300 fails: pkg = 80392 bytes (343x300 works)
