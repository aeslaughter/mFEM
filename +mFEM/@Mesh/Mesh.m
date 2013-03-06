classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access = protected)
        elements = Composite();%mFEM.elements.Line2.empty();
        nodes = Composite(); %mFEM.elements.base.Node.empty();
        n_nodes;                
        n_elements;
        options = struct('time', false, 'space', 'scalar');
    end
    
    properties (Access = protected)
        node_map = uint32([]);
        elem_map = uint32([]);
        boundary_map = logical([]);
        boundary_tag = {};
        dof_map = uint32([]);
    end
    
    methods
        plot(obj,varargin);
        grid(obj,varargin);
        
        function obj = Mesh(varargin)
            obj.options = gatherUserOptions(obj.options, varargin{:});  
        end
        
        function delete(obj)
            nodes = obj.nodes;
            elements = obj.elements;
            spmd
                nodes.delete();
                elements.delete();
            end
        end
        
        function node = createNode(obj, x)
            id = length(obj.nodes) + 1;
            node = mFEM.elements.base.Node(id,x);     
            obj.nodes{id} = node; 
        end
        
        function elem = createElement(obj, type, nodes)
            
            if isnumeric(nodes);
                nodes = obj.nodes(nodes);
            end
            
            id = length(obj.elements) + 1;
            elem = feval(['mFEM.elements.',type],id,nodes);
            obj.elements{id} = elem;
        end



        
        function elem = getElement(obj, id)
            elem = gather(obj.elements(id));
            elem = elem{1};
        end
    end
    
    methods (Access = protected)
        init(obj)        
        idEmptyBoundary(obj,id);

    end
    
    methods (Static)
        nodes = buildNodes(node_map);
        elements = buildElements(type,elem_map,node_map,nodes);
    end
end

