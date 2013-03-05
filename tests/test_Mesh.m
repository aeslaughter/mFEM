function test_Mesh

% Test 1D grid creation
mesh = mFEM.Mesh();
mesh.grid('Quad4',0,1,0,1,5,5);


% mFEM.elements.Quad4.grid(0,1,0,1,1,6)