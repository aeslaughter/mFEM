classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access = protected)
        elements = Composite();%mFEM.elements.Line2.empty();
        nodes = Composite(); %mFEM.elements.base.Node.empty();
        options = struct('time', false);
    end
    
    properties (Access = protected)
        node_map = uint32([]);
        elem_map = uint32([]);
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
        
        %addNode
        %addElement -> addNode
        
        function init(obj)

            obj.findNeighbors();
            % replace the numeric arrays with mFEM.Vectors
        end
        
        function elem = getElement(obj, id)
            elem = gather(obj.elements(id));
            elem = elem{1};
        end
        
%         function addElement(obj, elem_type, nodes)
%             % ADDELEMENT to FEmesh object (must be reinitialized)
%             %
%             % Syntax:
%             %   addElement(nodes);
%             %
%             % Input:
%             %   The input should take the exact form as the input for the
%             %   element that was established when the FEmesh object was
%             %   created. For example:
%             %
%             %   >> mesh = FEmesh();
%             %   >> mesh.addElement(1,'Quad4', [0,0; 1,0; 1,1; 0,1]);
%             %
%             %   Note, the optional space string for the element is not
%             %   accepted here, as it is defined across the entire mesh
%             %   when the FEmesh object is created.
%             
%             % Indicate that the object must be reinitialized
%             obj.initialized = false;
%                         
%             % Create the element
%             id = length(obj.element) + 1;
%             obj.element(id) = feval(['mFEM.elements.', elem_type],id, nodes);        
%         end
       
        

        
    end
    
    methods (Access = protected)
        findNeighbors(obj);
        

    end
    
    methods (Static)
        nodes = buildNodes(node_map);
        elements = buildElements(type,elem_map,node_map,nodes);
%         [node_map, elem_map] = buildGrid(order, varargin);
%         [node_map, elem_map] = buildMaps(varargin);
        
%         nodes = buildNodes(varargin);
%         [node_map, elem_map] = buildMaps(varargin)
        
    end
end

