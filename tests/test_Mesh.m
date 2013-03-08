function T = test_Mesh
%TEST_MESH Test function the the mFEM.Mesh class.

% Create an instance of mFEM.Test class
T = mFEM.Test();

% Generate a mesh
% try
    mesh = mFEM.Mesh('space','vector','-time');
    mesh.grid('Quad4',0,2,0,2,5,5);
%     mesh.plot([]);

% catch err
%     T.caught(err);
%     return;
% end

mesh.addBoundary('A','top','right');
mesh.addBoundary('B','x<1','y<1');
mesh.addSubdomain('C',{'x>1.5','y>1.5'});


mesh.getDof('Tag','A');
% mesh.getDof('Tag',{'A','B'});
    
    

% Test neighbor finding (test_Element has more extensive testing of this)
elem = mesh.getElements([3,15,23,11]);

% elem(1).sides(1).tag
% elem(1).nodes(1).tag
% elem(1).nodes(2).tag



delete(mesh);

