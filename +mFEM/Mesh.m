classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access = protected)
        elements = mFEM.elements.Line2.empty();
        nodes = mFEM.elements.base.Node.empty();
        options = struct('time', false);
    end
    
    methods
        function obj = Mesh(varargin)
            obj.options = gatherUserOptions(obj.options, varargin{:});  
        end
        
%         function delete(obj)
%            obj.elements.delete
%            obj.nodes.delete
%         end
        
        function node = createNode(obj, x)
            node = mFEM.elements.base.Node(x);     
            obj.nodes(end+1) = node; 
        end
        
        function elem = addElement(obj, type, nodes)
            id = length(obj.elements) + 1;
            elem = feval(['mFEM.elements.',type],id,nodes);
            obj.elements(id) = elem;
        end
        
        function init(obj)

            obj.findNeighbors();
            % replace the numeric arrays with mFEM.Vectors
        end
        
        function elem = getElement(obj, id)
            elem = gather(obj.elements.get(id));
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
        
        function grid(obj, type, varargin)
            
            % Display wait message
            if obj.options.time;
                ticID = tMessage('Generating Grid...');
            end
            
            [obj.nodes, obj.elements] = ...
                feval(['mFEM.elements.',type,'.grid'], varargin{:});

            % Complete message
            if obj.options.time;
                tMessage(ticID);
            end; 
        end
    end
    
    methods (Access = protected)
        function findNeighbors(obj)
            % Locates elements that share a side
            %
            % Syntax
            %   findNeighbors()   
            %
            % Description
            %   findNeighbors() locates neighboring elements and
            %   neighboring sides for the finite element mesh, it sets the
            %   various neighbor and side related properties for the
            %   elements in the mesh.

            % Display wait message
            if obj.options.time;
                ticID = tMessage('Locating neighbor elements...');
            end
            
            % Loop through each element and find neighboring elements
            elements = obj.elements;
            
            spmd
               local = getLocalPart(elements);
               
               for e = 1:length(local)
                  elem = local{e};
                  
                  neighbors = elem.getNodeParents();
                  
                  for i = 1:length(neighbors);
                      for s = 1:length(elem.sides);
                          
                          if ~isempty(elem.sides(s).neighbor); 
                              break;
                          end
                          
                      	  [C,~,ns] = intersect(elem.side_map(s,:), neighbors(i).side_map, 'rows');
                          if ~isempty(C);
                             elem.sides(s).neighbor = neighbors(i);
                             elem.sides(s).neighbor_side = ns;
                             neighbors(i).sides(ns).neighbor = elem;
                             neighbors(i).sides(ns).neighbor_side = s;   
                          else
                             elem.sides(s).neighbor = [];
                             elem.sides(s).neighbor_side = [];
                          end
                      end
                  end
               end 
            end
                      
            % Complete message
            if obj.options.time;
                tMessage(ticID);
            end;
        end 
    end
    
%     methods (Static)
%     end
end

