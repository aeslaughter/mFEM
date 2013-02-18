classdef Cell < handle & matlab.mixin.Heterogeneous
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties %(GetAccess = public, SetAccess = ?mFEM.Mesh)
        id = [];
        n_nodes = [];
        nodes = {};
        sides = struct('neighbor',{},'neighbor_side',{});      
    end   
    
    properties (Abstract, Access = protected) 
       side_ids;           
       n_sides;
    end
    
    properties %(Access = ?mFEM.Mesh)
        side_map = mFEM.elements.base.Node.empty();%
        node_parent_map = mFEM.elements.Line2.empty();%.
    end
    
    methods (Abstract)
       % build(obj);
       %refine(obj);
       %coarsen(obj);
       %plot(obj)
    end
    
    methods
        function obj = Cell(id, nodes, varargin)
           obj.nodes = nodes;
           obj.id = id;
           
           for i = 1:obj.n_sides;
               obj.side_map(i,:) = obj.nodes{obj.side_ids(i,:)};
           end

           for i = 1:length(nodes);
               nodes{i}.addParent(obj);
           end
           
           obj.sides(obj.n_sides) = struct('neighbor',[],'neighbor_side',[]);
        end
        
        function out = getNodeParents(obj)
            
           out = mFEM.elements.Line2.empty();
           for i = 1:length(obj.nodes);
                out = [out, obj.nodes{i}.parents(obj.nodes{i}.parents ~= obj)];
           end
           out = unique(out);
           
        end
    end
    
end

