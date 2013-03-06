classdef Node < handle
    %Node Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coord = [0;0;0];
        id = [];
        on_boundary = false;
        n_dim;
        dof = uint32([]);
        tag = {};
    end
    
    properties (GetAccess = public, SetAccess = protected)
        parents;
    end
    
    methods
        function obj = Node(id,x)
            if nargin == 2;
            	obj.init(id,x);
            end
        end
        
        function init(obj,id,x)           
            if iscell(x);
                for i = 1:length(obj);
                    n = length(x{i});
                    obj(i).id = id(i);
                    obj(i).n_dim = n;
                    obj(i).coord(1:n,1) = x{i};
                end   
            else
                n = size(x,2);
                for i = 1:length(obj);
                    obj(i).id = id(i);
                    obj(i).n_dim = n;
                    obj(i).coord(1:n,1) = x(i,:);
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

