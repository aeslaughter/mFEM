function test_Mesh

% profile on;
mesh = mFEM.Mesh('-time','space','vector');
mesh.grid('Quad4',0,2,0,2,3,3);

mesh.addTag('A','x<=1');
mesh.addTag('B','right');


delete(mesh);

% 344x300 fails: pkg = 80392 bytes (343x300 works)
