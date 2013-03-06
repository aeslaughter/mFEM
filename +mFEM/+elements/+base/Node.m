classdef Node < handle
    %Node Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coord = [0;0;0];
        id = [];
        on_boundary = true;
    end
    
    properties (GetAccess = public, SetAccess = protected)
        parents;
    end
    
    methods
        function obj = Node(id, x)
            if nargin == 2;
            	obj.init(id,x');
            end
        end
        
        function init(obj,id,x)
            n = size(x,2);
            for i = 1:length(obj);
                obj(i).id = id(i);
                obj(i).coord(1:n,1) = x(i,:);
            end
        end
        
        function out = getCoord(obj)
            out = [obj.coord]';
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
        
        
    end
    
    methods (Access = ?mFEM.elements.base.Element)
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
    end
end

