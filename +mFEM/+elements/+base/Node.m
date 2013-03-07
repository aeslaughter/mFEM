classdef Node < handle
    %NODE Class for defining node objects
    % Inludes the general behavior for a finite element node, including the 
    % node, global degrees-of-freedom.   
    %
    % Syntax
    %   nodes = mFEM.elements.base.Node(id,coord)
    %   nodes(n) = mFEM.elements.base.Node()
    %
    % 
    % Description (todo)
    %   nodes = mFEM.elements.base.Node(id,coord)
    %   nodes(n) = mFEM.elements.base.Node()
    %
    % See Also mFEM.elements.base.Element
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
        n_dim;                  % no. of spatial dimensions
        dof = uint32([]);       % global degrees-of-freedom
        tag = {};               % list of char tags for this node
        lab = uint32([]);       % the processor that holds this node
        parents;                % handles of elements with this node
    end
    
    % Begin method definitions
    methods
        function obj = Node(id,x)
            if nargin == 2;
            	obj.init(id,x);
            end
        end
        
        function init(obj,id,x,varargin) 
            
           lab = 1; 
           if nargin == 4 && isinteger(varargin{1});
               lab = varargin{1};
           end
           
            if iscell(x);
                for i = 1:length(obj);
                    n = length(x{i});
                    obj(i).id = id(i);
                    obj(i).n_dim = n;
                    obj(i).coord(1:n,1) = x{i};
                    obj(i).lab;
                end   
            else
                n = size(x,2);
                for i = 1:length(obj);
                    obj(i).id = id(i);
                    obj(i).n_dim = n;
                    obj(i).coord(1:n,1) = x(i,:);
                    obj(i).lab;
                end
            end
        end
        
        function out = getCoord(obj)
            out = [obj.coord];
        end
        
        function addTag(obj,tag)
           for i = 1:length(obj);
               n = length(obj(i).tag)+1;
               obj(i).tag{n} = tag;
           end
        end
        
        function tf = hasTag(obj,tag)
            n = length(obj);
            tf = false(n,1);
            for i = 1:n;
                tf(i) = any(strcmp(tag,obj(i).tag));
            end
        end
        
        function varargout = get(obj, varargin)
           
            if nargin == 1
                if nargout == 1
                    varargout{1} = obj.coord;
                else
                    varargout = num2cell(obj.coord);
                end
            
            elseif isnumeric(varargin{1});
                varargout{1} = obj.coord(varargin{1});
                
            elseif ischar(varargin{1});
                switch lower(varargin{1});
                    case 'x'; varargout{1} = obj.coord(1);
                    case 'y'; varargout{1} = obj.coord(2);
                    case 'z'; varargout{1} = obj.coord(3);
                end
            end  
        end
        
%         function getDof(obj,varargin)
%             for i = 1:length(obj);  
%             end
%         end
    end
    
    methods (Access = {?mFEM.elements.base.Element,?mFEM.Mesh})
        function addParent(obj, elem)
            for i = 1:length(obj);
                idx = length(obj(i).parents);
                if idx == 0;
                    obj(i).parents = elem;
                else
                    obj(i).parents(idx+1) = elem;
                end
            end
        end 

        function out = getParents(obj,varargin)
            out = {};
            for i = 1:length(obj);
                out = [out,obj(i).parents];
            end
            out = unique(out);
            
            if nargin == 2;
                out = out(out~=varargin{1});
            end
        end
        
        function setDof(obj,type)
            for i = 1:length(obj);
                if strcmpi(type,'vector'); 
                    n = obj(i).n_dim; 
                else
                    n = 1;
                end
                obj(i).dof = transformDof(obj(i).id,n);
            end
        end
        
        function setBoundaryFlag(obj)
           for i = 1:length(obj);
               obj(i).on_boundary = true;
           end
        end
    end
end

