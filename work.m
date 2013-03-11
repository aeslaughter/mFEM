function work

mesh = mFEM.Mesh();
mesh.grid('Line2',0,4,2); 

elem = mesh.getElements();
no = mesh.getNodes('-gather');



