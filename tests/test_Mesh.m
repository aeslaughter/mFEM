function T = test_Mesh
    %TEST_MESH Test function the the mFEM.Mesh class.
    %
    %     13---14----15---16
    %     |  7  |  8  |  9 |
    %     9----10----11---12
    %     |  4  |  5  |  6 | 
    %     5-----6-----7----8
    %     |  1  |  2  |  3 | 
    %     1-----2-----3----4
    
    % Create an instance of mFEM.Test class
    T = mFEM.Test();

    % Generate a mesh
    % try
        mesh = mFEM.Mesh('space','scalar','-time');
        mesh.grid('Quad4',0,3,0,3,3,3);
    %     mesh.plot([]);

    % catch err
    %     T.caught(err);
    %     return;
    % end
% 
%     mesh.addBoundary('A','top','right');    
%     mesh.addBoundary('B','x<1','y<1');
%     mesh.addSubdomain('C',{'x>1.5','y>1.5'});
%     % add.init()
% 
%     no = mesh.getNodes();
%     el = mesh.getElements();
%     
% 
%     no(4)
%     el(3).nodes(2)




% dof1 = mesh.getDof('Tag','A','-index');
% dof2 = mesh.getDof('Tag','A','-gather');
% 
% idx(:,1) = 1:mesh.n_dof;
% idx(dof2) == dof1;

% mesh.getDof('Tag',{'A','B'});
    
    

% Test neighbor finding (test_Element has more extensive testing of this)

% elem(1).sides(1).tag
% elem(1).nodes(1).tag
% elem(1).nodes(2).tag

delete(mesh);

