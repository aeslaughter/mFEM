function test_Mesh

profile on;
mesh = mFEM.Mesh('-time');
mesh.grid('Quad4',0,2,0,2,20,20);
% mesh.plot([]);
profile viewer;


delete(mesh);

% 344x300 fails: pkg = 80392 bytes (343x300 works)
