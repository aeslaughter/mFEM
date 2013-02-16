classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access = protected)
        elements = mFEM.elements.Line2.empty();
        nodes = mFEM.elements.base.Node.empty();
        map = struct('node',{});
    end
    
    methods
        function obj = Mesh()
            
        end
        
        function node = addNode(obj, x)
            node = mFEM.elements.base.Node(x);     
            obj.nodes(end+1) = node; 
        end
        
        function elem = addElement(obj, type, nodes)
            id = length(obj.elements) + 1;
            elem = feval(['mFEM.elements.',type],id,nodes);
            obj.elements(id) = elem;
        end
        
        function init(obj)
            % distribute nodes and elements (if not already, grid) and
            % replace the numeric arrays with mFEM.Vectors
        end
        
%         function elem = getElement(obj, id)
%             vec = obj.elements;
%             
%             spmd  
%                 local = getLocalPart(vec);
%                 idx = [];
%                 for i = 1:length(local);
%                     if local{i}.id == id;
%                         idx = i;
%                         break;
%                     end
%                 end 
% 
%                 codist = codistributor1d(codistributor1d.unsetDimension, ...
%                     codistributor1d.unsetPartition, [numlabs,1]);
%                 IDX = codistributed.build(idx, codist); 
%             end
%             
%             c = gather(vec(IDX));
%             elem = c{1};
%         end
        
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
        
        function grid(obj, cell_type, varargin)
            
            % cell_type = LINE, TRI, QUAD, TET, HEX, PRISIM
            
            n = length(varargin);
            if n == 3 || (n > 3 && ischar(varargin{4}))
                obj.gen1Dgrid(elem_type, varargin{:});

%             elseif n == 6 || (n > 6 && ischar(varargin{7}))
%                 obj.gen2Dgrid(elem_type, varargin{:});
% 
%             elseif n == 9 || (n > 9 && ischar(varargin{10}))
%                 obj.gen3Dgrid(elem_type, varargin{:});

            else
                error('Mesh:grid:InputError', 'The type of grid desired was not understood, check that the proper number of inputs were supplied.');
            end  
        end
    end
    
    methods (Access = protected)
                  
        function gen1Dgrid(obj, elem_type, x0, x1, xn, varargin)

            dx = (x1-x0)/xn;
            nX = x0 : dx : x1;
            eX = 1:length(nX)-1;
            node_global_size = [length(nX),1];
            elem_global_size = [length(eX),1];
            nX = distributed(nX);
            eX = distributed(eX);

            spmd 
                x = getLocalPart(nX);
                id = globalIndices(eX,1);
                
                node_local = cell(length(x),1);
                elem_local = cell(length(id),1);
                for i = 1:length(x);
                    node_local{i} = mFEM.elements.base.Node(x(i));
                    
                    if i > 1;
                        elem_local{i-1} = feval(['mFEM.elements.', elem_type], id(i-1), [local_node{i-1},local_node{i}]);
                    end
                end
    
                 node_codist = codistributor1d(codistributor1d.unsetDimension, ...
                    codistributor1d.unsetPartition, node_global_size);
                 elem_codist = codistributor1d(codistributor1d.unsetDimension, ...
                    codistributor1d.unsetPartition, elem_global_size);
                 node_dist = codistributed.build(node_local, node_codist);
                 elem_dist = codistributed.build(elem_local, elem_codist);
            end
            obj.nodes = mFEM.Vector(node_dist);
            obj.elements = mFEM.Vector(elem_dist);
        end
    end
end

