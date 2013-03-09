function T = test_Node(varargin)
    %TEST_NODE Tests all methods of the node class.
    %   This test creates 6 nodes with random 2D coordinates and evaluates
    %   and test the output of all method, public and protected of the Node
    %   class.
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
    
    % Create the test class
    T = mFEM.Test('name','Node',varargin{:});   
    
    % Create the node objects (use empty construtor method)
    try
        n = 6;   
        coord = rand(n,2);
        node(n) = mFEM.elements.base.Node();
    catch err
        T.caught(err);
        return;
    end

    % Intialize the nodes
    T.handle = node; % use test class to perform methods
    T.eval('init',1:n,coord,3); % initialize
    T.compare(node(4).getCoord(),coord(4,:)','Extracted coordinates');
    
    % Test dof set/get
    T.eval('setDof',1); % set 3 dofs per node
    T.handle = node(6);
    dof = T.eval('getDof');
    T.compare(dof,[16,17,18],'Extracted completed dofs from node 6');
    dof = T.eval('getDof','Component',[2,3]);
    T.compare(dof,[17,18],'Extracted 2,3 dof component from node 6');
    dof = T.eval('getDof','Component','x');
    T.compare(dof,16,'Extracted x-component dof from node 6');
    
    % Set the boundary id for nodes 3 & 6
    T.handle = node([3,6]);
    T.eval('setBoundaryFlag');
    T.compare(node(3).on_boundary,true,'Node 3 on boundary');
    T.compare(node(6).on_boundary,true,'Node 6 on boundary');
    
    % Add/get parent to nodes 5 and 2
    T.handle = node([2,5]);
    T.eval('addParent',2);
    p = T.eval('getParents');
    T.compare(p,2,'Parents extracted from nodes');

    % Test the tag
    T.eval('addTag','A');
    T.compare(node(5).tag{1},'A','Tag added correctly');
    
    % Complete the test
    delete(T);