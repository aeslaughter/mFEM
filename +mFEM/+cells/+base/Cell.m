classdef Cell < handle & matlab.mixin.Heterogeneous
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        level = 0;
        order = 1;
        children;
        vertices;
    end
    
    properties (Abstract)
       %side; 
       n_dim;
    end
    
    properties (Abstract, Access = protected)
        has_order;
    end
    
    methods (Abstract)
       % build(obj);
       %refine(obj);
       %coarsen(obj);
       %plot(obj)
    end
    
    methods
        function obj = Cell(vertices, varargin)
           obj.vertices = vertices; 
           
           if nargin == 2;
               if all(obj.has_order ~= varargin{1});
                   error('Cell:Cell:InvalidOrder', 'The supplied order of %d does not exist for this cell type.',varargin{1});
               end
               obj.order = varargin{1};
           end
        end
        

    end
    
end

