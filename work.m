function work

mesh = mFEM.Mesh();
mesh.grid('Line2',0,1,2);
% mesh.addBoundary('test','right');
% mesh.update();
elem = mesh.getElements(1,'-gather');

elem.shape(0)
elem.shapeDeriv(0)

% delete(mesh)
% delete(elem)
clear mesh elem clear
