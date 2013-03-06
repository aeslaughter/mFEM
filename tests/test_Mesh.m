function test_Mesh

% profile on;
mesh = mFEM.Mesh('-time','space','vector');
mesh.grid('Quad4',0,2,0,2,20,10);
% mesh.plot([]);
% profile viewer;

n1 = mesh.nodes{1}; n1(132).dof
delete(mesh);

% 344x300 fails: pkg = 80392 bytes (343x300 works)
