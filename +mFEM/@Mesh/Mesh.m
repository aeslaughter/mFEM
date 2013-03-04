classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access = protected)
        elements = mFEM.elements.Line2.empty();
        nodes = {}; %mFEM.elements.base.Node.empty();
        options = struct('time', false);
    end
    
    properties (Access = protected)
        node_map = uint32([]);
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
        
        function grid(obj, type, varargin)
            
            % Display wait message
            if obj.options.time;
                ticID = tMessage('Generating Grid...');
            end

            profile on;
            [node_map, elem_map] = feval(['mFEM.elements.',type,'.buildMaps'],varargin{:});

            nodes = obj.buildNodes(node_map);
            
            
            
            
            
%             feval(['mFEM.elements.',type,'.grid'],varargin{:});

            % Complete message
            if obj.options.time;
                tMessage(ticID);
            end; 
        end
        
        function plot(obj, data, varargin)
            
            opt.data = data;
            opt.labelelements = false;
            opt.labelnodes = false;
            opt = gatherUserOptions(opt,varargin{:});

            p = zeros(length(obj.elements),1);

            figure; hold on;
            elements = gather(obj.elements);
            for e = 1:length(elements);
                elem = elements{e};
                no = elem.getNodeCoord();
                face = elem.side_ids;
                p(e) = elem.plot(no,face);
                
                if opt.labelelements;
                    cntr = num2cell(mean(no,1));
                    text(cntr{:}, num2str(elem.id),'FontSize',14,...
                    'BackgroundColor','k','FontWeight','Bold',...
                    'Color','w','HorizontalAlignment','center');
                end
                
                if opt.labelnodes;
                    for i = 1:length(elem.nodes);
                        n = elem.nodes{i};
                        loc = num2cell(n.coord);
                        text(loc{:},num2str(n.id),'FontSize',10,'Color','w',...
                            'BackgroundColor','b','HorizontalAlignment','center');
                    end
                end
            end
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
                      for s = 1:length(elem.side_ids);
                          
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
    
    methods (Static)
        nodes = buildNodes(varargin);
%         [node_map, elem_map] = buildGrid(order, varargin);
%         [node_map, elem_map] = buildMaps(varargin);
        
%         nodes = buildNodes(varargin);
%         [node_map, elem_map] = buildMaps(varargin)
%         elements = buildElements(nodes, elem_map);
    end
end

