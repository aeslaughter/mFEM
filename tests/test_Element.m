function T = test_Element(varargin)
    %TEST_ELEMENT Tests the general behavior of an element
    %   The following is the mesh that is used for the testing, which runs
    %   the non-shape function related methods of the Element class, using
    %   the Quad4 element.
    %
    %     5-----6-----7----8
    %     |  1  |  2  |  3 | 
    %     1-----2-----3----4
    %
    % Lab 1: Nodes = [1,2,3,4], Elements = [1,2]
    % Lab 2: Nodes = [5,6,7,8], Elements = [3]
    %
    % See Also bin\test mFEM.Test
    %
    %----------------------------------------------------------------------
    %  mFEM: A Parallel, Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2013 Andrew E Slaughter
    % 
    %  This program is free software: you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation, either version 3 of the License, or
    %  (at your option) any later version.
    % 
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    % 
    %  You should have received a copy of the GNU General Public License
    %  along with this program. If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
   
    % Create the test object
    T = mFEM.Test('Name','Element',varargin{:});
    
    % Start a parallel operation
    openpool = false;
    if matlabpool('size') == 0;
        matlabpool(2);
        openpool = true;
    elseif matlabpool('size') ~= 2;
        T.result(false,'Test only operates in parallel, with two labs');
    end
    
    % Create the mesh
    try
        mesh = mFEM.Mesh('Space','Vector');
        mesh.grid('Quad4',0,3,0,3,3,1);  
        mesh.addBoundary(1,'top');
        mesh.update();
    catch err
        T.caught(err);
        return
    end
    
    % Extract the nodes and elements from the two labs
    try
        N1 = mesh.getNodes('lab',1);
        N2 = mesh.getNodes('lab',2);
        E1 = mesh.getElements('lab',1);
        E2 = mesh.getElements('lab',2);
    catch err
        T.caught(err)
        return
    end

    % Pefrom direct node vs. element accessed node dof tests
    T.compare(N1(3).dof, E1(2).nodes(2).dof,'Degrees-of-freedom extraction, lab 1');
    T.compare(N1(3).dof, E2(1).nodes(1).dof,'Degrees-of-freedom extraction, lab 2');
    T.compare(N2(3).getDof('Component',2), E1(2).nodes(3).getDof('Component',2),'Degrees-of-freedom access, limited to y-component, lab 1');
    T.compare(N2(3).getDof('Component',2), E2(1).nodes(4).getDof('Component',2),'Degrees-of-freedom access, limited to y-component, lab 2');

    % Test tags
    T.compare(N2(3).tag{1}, E1(2).nodes(3).tag{1},'addBoundary tag access, lab 1');
    T.compare(N2(3).tag{1}, E2(1).nodes(4).tag{1},'addBoundary tag access, lab 2');
 
    % Test side dofs (across labs)
    T.compare(sort(E1(2).getDof('side',2)),sort(E2(1).getDof('side',4)),'Element side dofs captured');
    T.compare(E1(1).getDof('side',3),[11,12,9,10],'Element side global dofs accsses correctly');
    T.compare(E1(1).getDof('side',3,'-local'),[5,6,7,8],'Element side local dofs accessed correctly');
    
    % Clean up after the testing
    clear mesh N1 N2 E1 E2 T
    if openpool;
        matlabpool close;
    end
end