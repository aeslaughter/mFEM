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
    properties %(GetAccess = public, SetAccess = ?mFEM.Mesh)
        coord = [0;0;0];        % column vector of spatial coordinates
        id = uint32([]);        % unique global id
        on_boundary = false;    % true if node is on a boundary
        n_dim;                  % no. of spatial degrees-of-freedom
        n_dof;                  % no. of degrees-of-freedom for node
        dof = uint32([]);       % global degrees-of-freedom
        tag = {};               % list of char tags for this node
        lab = uint32([]);       % the processor that holds this node
        parents;                % handles of elements with this node
    end
    
    % Begin method definitions
    methods
        function obj = Node(varargin)
            if nargin > 0;
            	obj.init(varargin{:});
            end
        end
        
        function init(obj,id,x,varargin) 
            %INIT Initializes the node objects 
            
            ndim = size(x,2);
            % numeric | 'scalar' | 'vector'
           if nargin == 4 && isnumeric(varargin{1});
               ndof = varargin{1};
           elseif nargin == 4 && ischar(varargin{1}) && strcmpi('vector',varargin{1});
               ndof = ndim;
           else
               ndof = 1;
           end
             
 
           n = length(obj);
           x = mat2cell(x',ndim,ones(1,n));
           ndim = num2cell(repmat(ndim,n,1));
           ndof = num2cell(repmat(ndof,n,1));
           proc = num2cell(repmat(labindex,n,1));
           id = num2cell(id);
       
           [obj.id] = id{:};
           [obj.n_dim] = ndim{:};
           [obj.n_dof] = ndof{:};
           [obj.lab] = proc{:};
           [obj.coord] = x{:};

% Loop was about slower (29s vs. 18s for 1 mil. nodes)
%             lab = labindex;
%             for i = 1:length(obj);
%                 obj(i).id = id(i);
%                 obj(i).n_dim = n_dim;
%                 obj(i).n_dof = n_dof;
%                 obj(i).coord(1:n_dim,1) = x(i,:);
%                 obj(i).lab = lab;
%             end
        end
        
        function out = getCoord(obj)
            out = [obj.coord];
        end
        
        function dof = getDof(obj,varargin)
            
            opt.component = [];
            opt = gatherUserOptions(opt,varargin{:});
            cmp = obj.parseComponentInput(opt.component);
            
             if isempty(cmp);
                 dof = [obj.dof];  
             else
                dof = {obj.dof};
                 for i = 1:length(dof);
                     dof{i} = dof{i}(cmp);
                 end
                 dof = cell2mat(dof); 
             end
        end
        
        function addTag(obj,tag)
           for i = 1:length(obj);
               n = length(obj(i).tag)+1;
               obj(i).tag{n} = tag;
           end
        end
        
%         function tf = hasTag(obj,tag)
%             n = length(obj);
%             tf = false(n,1);
%             for i = 1:n;
%                 tf(i) = any(strcmp(tag,obj(i).tag));
%             end
%         end
        
%         function varargout = get(obj, varargin)
%            
%             if nargin == 1
%                 if nargout == 1
%                     varargout{1} = obj.coord;
%                 else
%                     varargout = num2cell(obj.coord);
%                 end
%             
%             elseif isnumeric(varargin{1});
%                 varargout{1} = obj.coord(varargin{1});
%                 
%             elseif ischar(varargin{1});
%                 switch lower(varargin{1});
%                     case 'x'; varargout{1} = obj.coord(1);
%                     case 'y'; varargout{1} = obj.coord(2);
%                     case 'z'; varargout{1} = obj.coord(3);
%                 end
%             end  
%         end
        
%         function getDof(obj,varargin)
%             for i = 1:length(obj);  
%             end
%         end
    end
    
    methods %(Access = {?mFEM.elements.base.Element,?mFEM.Mesh)
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
        
        function dof = setDof(obj,varargin)
            %SETDOF sets the degrees-of-freedom for the node(s)
            %
            % Syntax
            %   setDof()
            %   setDof(strt)
            %   setDof(dof)
            %
            % Description
            %
            % dof(n,:), n = size(obj)  (implicit)
            % dof = scalar, strt index (explicit)
            
            if nargin == 2;
                dof = varargin{1};
            else
                dof = 1;
            end
            
            if isscalar(dof);
                strt = dof;
                for i = 1:length(obj);
                    stop = strt + obj(i).n_dof - 1;
                    obj(i).dof = strt:stop;
                    obj(i).dof
                    strt = stop+1;
                end
            elseif size(dof,1) == length(obj);
                for i = 1:length(obj);
                    obj(i).dof = dof(i,:);
                end
            end
        end

        function setBoundaryFlag(obj)
           for i = 1:length(obj);
               obj(i).on_boundary = true;
           end
        end
    end
    
    methods (Static)
        function cmp = parseComponentInput(cmp)
            
            msgid = 'Node:getDof:InvalidComponent';
            errmsg = 'The supplied component must be number or ''x'', ''y'', or ''z''.';
            
            if ischar(cmp); cmp = {cmp}; end
             if iscell(cmp);
                 for i = 1:length(cmp);
                     if ischar(cmp{i});
                         switch lower(cmp{i});
                             case 'x'; cmp{i} = 1;
                             case 'y'; cmp{i} = 2;
                             case 'z'; cmp{i} = 3;
                             otherwise
                                error(msgid,errmsg);
                         end
                     elseif ~isnumeric(cmp{i});
                        error(msgid,errmsg);
                     end
                 end
                 cmp = cell2mat(cmp);
             elseif ~isnumeric(cmp);
                 error(msgid,errmsg);
             end
        end
    end
end

