classdef Element < mFEM.elements.base.HideHandle
    %ELEMENT Base class for defining elements.
    %   Inludes the general behavior of an element, including the node 
    %   locations, id, shape functions, etc...
    %
    %   This is an abstract class, so it must be inherited to function.
    %   The abstract properties and methods must be redifined in the
    %   subclass, see Line2.m for an example. In general, if you need help 
    %   for an element see the the help for the subclass itself.
    %
    %   IMPORTANT: The nodes stored by the element are copies of the nodes
    %   stored in the mFEM.Mesh property. If node operations are performed
    %   (i.e., addBoundary or addSubdomain) the update() method of the Mesh
    %   class must be called for those changes to be reflected on the
    %   element copies. This only affects parallel operations where the
    %   elements contain nodes from other processors.
    %
    % See Also Quad4
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
    properties %(GetAccess = public, SetAccess = {?mFEM.Mesh,?mFEM.elements.base.Element})
        id = uint32([]);        % unique global id
        on_boundary = false;    % true when element touches a border
        sides;                  % structure of side information
        tag = {};               % list of char tags for this element
        n_dof;                  % total no. of dofs associated with element
        dof;                    % dofs for the element, in order of nodes
        n_sides;                % no. of sides for this element
        n_nodes;                % no. of nodes for this element
        nodes =...              % node objects for this element
            mFEM.elements.base.Node.empty();
    end
    
    % Constants that must be defined by inhering class (e.g., Quad4)
    properties (Abstract, Constant, Access = public)
        side_ids;   % map that defines which nodes comprise each side  
        side_type;  % name of element for sides
        quad;       % Gauss quadrature rules to utilize
        n_dim;      % no. of spacial dimensions
    end

    properties (Access = ?mFEM.Mesh)
        lab;                    % the processor that holds this element
    end
    
    % Abstract Methods (protected)
    % (the user must redfine these in subclasse, e.g. Line2)
    methods (Abstract, Access = protected)
        N = basis(obj, varargin)            % basis functions
        B = gradBasis(obj, varargin)        % basis function derivatives (dN/dx, ...)
        G = localGradBasis(obj, varargin)   % basis function derivatives (dN/dxi, ...)
        J = jacobian(obj, varargin)         % the Jacobian matrix for the element
    end

    % Public Methods
    % (These methods are accessible by the user to create the element and
    % access the shape functions and other necessary parameters)
    methods (Access = public)
        init(obj,id,nodes);
        dof = getDof(obj,varargin);
        nodes = getNodes(obj);
        varargout = hasTag(obj,tag);
        N = shape(obj,x,varargin);
        B = shapeDeriv(obj,x);
        J = detJ(obj,x);
%         delete(obj);

        function obj = Element(varargin)
            %ELEMENT Class constructor.
            %
            % This is an abstract class, it must be inherited by a subclass
            % to operate, see Line2.m for example. The following syntax and
            % descriptions apply to all subclasses unless noted otherwise
            % in the documentation for the specific element.
            %
            % Syntax
            %   Element(id, nodes)
            %   Element(id, nodes, 'PropertyName', PropertyValue, ...)
            %
            % Description
            %   Element(id, nodes) creates an element given, where id is a
            %   unique identification number for this element and nodes is 
            %   a matrix of node coordinates (global) that should be 
            %   arranged as column matrix (no. nodes x no. dims).
            %
            %   Element(id, nodes, 'PropertyName', PropertyValue, ...) 
            %   allows the user to customize the behavior of the element, 
            %   the available properties are listed below.

            if nargin > 0;
                obj.init(varargin{:});
            end
        end
    end
    
    methods %(Access = ?mFEM.Mesh)
         findNeighbors(obj);
         addTag(obj,tag,type);
         setDof(obj);
    end

    methods (Static)
        function buildNodeMap(varargin)
            error('Element:buildNodeMap:NotImplemented', 'The ''buildNodeMap'' method is not defined for this element, add the method to the parent class (e.g., Quad4.m)');
        end
        function buildNodess(varargin)
            error('Element:buildElementMap:NotImplemented', 'The ''buildElementMap'' method is not defined for this element, add the method to the parent class (e.g., Quad4.m)');
        end
    end
end