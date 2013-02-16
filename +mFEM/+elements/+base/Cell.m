classdef Cell < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dim;
        level = 0;
        order = 1;
        children = mFEM.elements.base.Cell.empty();
        nodes = mFEM.elements.base.Node.empty();
    end
    
    properties (Abstract)
       side; 
    end
    
    methods
        function obj = Cell(dim,nodes,varargin)
           obj.dim = dim;
           obj.nodes = nodes; 
           
           if nargin == 3;
               obj.level = varargin{1};
           end
        end
        
        function refine(obj)
            
            x = obj.getNodeCoord();
            
            n = mFEM.elements.base.Node(mean(x));
            
            obj.children(1) = mFEM.elements.base.Cell(obj.dim,[obj.nodes(1),n],obj.level+1);
            obj.children(2) = mFEM.elements.base.Cell(obj.dim, [n,obj.nodes(2)],obj.level+1);
            
            
            
        end
        
        function [x,y,z] = getNodeCoord(obj)
           
            for i = 1:length(obj.nodes);
               [x(i,:),y(i,:),z(i,:)] = obj.nodes(i).get(); 
            end
            
            
        end
        
    end
    
end

