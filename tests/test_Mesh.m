function test_Mesh

mesh = mFEM.Mesh();
tic;
mesh.grid('Quad4',0,1,0,1,200,200); 
toc;


% 344x300 fails: pkg = 80392 bytes (343x300 works)
