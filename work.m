function work

mesh = mFEM.Mesh();
mesh.grid('Line2',0,1,2);
elem = mesh.getElements(1,'-gather');
mesh.addBoundary('test','right');

elem.shape(0)
elem.shapeDeriv(0)

delete(mesh)
delete(elem)
clear mesh elem clear
