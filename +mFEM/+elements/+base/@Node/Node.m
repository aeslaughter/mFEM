classdef Node < mFEM.elements.base.HideHandle
    %NODE Class for defining node objects
    %   Inludes the general behavior for a finite element node, including 
    %   coordinates, boundary information, and global degrees-of-freedom.
    %
    %   The class is setup to allow for the handling of nodes with varying
    %   spatial dimensions and no. of degrees-of-freedom per node.
    %
    % Syntax
    %   nodes = mFEM.elements.base.Node(id,coord)
    %   nodes(n) = mFEM.elements.base.Node()
    %
    % Description
    %   nodes = mFEM.elements.base.Node(id,coord) this is the standard
    %   method for creating nodes individually, which is slow. The id is
    %   the global identifier (integer) and the coordinates should be a
    %   column vector of nodal position.
    %
    %   nodes(n) = mFEM.elements.base.Node() used the empty constructor to
    %   create an array of elements and then call the init method to setup
    %   the nodes. This method is significantly faster (order of magnitude)
    %   than looping with the first method.
    %
    % See Also init Element
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
    
    % Basic read-only properties that are available for user
    properties (GetAccess = public, SetAccess = ?mFEM.Mesh)
        coord = [0;0;0];        % column vector of spatial coordinates
        id = uint32([]);        % unique global id
        on_boundary = false;    % true if node is on a boundary
        n_dim;                  % no. of spatial degrees-of-freedom
        n_dof;                  % no. of degrees-of-freedom for node
        dof = uint32([]);       % global degrees-of-freedom
        tag = {};               % list of char tags for this node
        lab = uint32([]);       % the processor that holds this node
    end
    
    properties (Access = {?mFEM.elements.base.Element,?mFEM.Mesh})
        parents;                % handles elements with this node
    end
    
    % Public method definitions
    methods
        init(obj,id,x,varargin); 
        out = getCoord(obj);
        dof = getDof(obj,varargin);
        addTag(obj,tag);
        
        function obj = Node(varargin)
            % Class constructor (see help mFEM.elements.base.Node)
            if nargin > 0;
            	obj.init(varargin{:});
            end
        end
    end
    
    % Protected methods
    methods (Access = {?mFEM.elements.base.Element,?mFEM.Mesh,?mFEM.Test})
        setDof(obj,varargin);
        setBoundaryFlag(obj);
        addParent(obj,elem);
        elem = getParents(obj,varargin);       
%         function resetParents(obj)
%             for i = 1:length(obj);
%                 obj(i).parents = [];
%             end
%         end
    end
   
    % Static methods
    methods (Hidden, Static, Access = protected)
        cmp = parseComponentInput(cmp);
    end
end