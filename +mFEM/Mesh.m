classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access = protected)
        elements
        
        map = struct('node',{});
    end
    
    methods
        function obj = Mesh()
            
        end
        
        function elem = getElement(obj, id)
            vec = obj.elements;
            
            spmd  
                local = getLocalPart(vec);
                idx = [];
                for i = 1:length(local);
                    if local{i}.id == id;
                        idx = i;
                        break;
                    end
                end 

                codist = codistributor1d(codistributor1d.unsetDimension, ...
                    codistributor1d.unsetPartition, [numlabs,1]);
                IDX = codistributed.build(idx, codist); 
            end
            
            c = gather(vec(IDX));
            elem = c{1};
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
        
        function grid(obj, elem_type, varargin)
            
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
        % Collect user options
%             opt.pol2cart = false;
%             opt.element = obj.options.element;
%             opt = gatherUserOptions(opt,varargin{:});
            
            % Generate the generic grid points
            dx = (x1-x0)/xn;
            X(1,:) = x0 : dx : x1-dx;
            X(2,:) = X(:,1) + dx;
            global_size = [length(X),1];
            X = distributed(X);
            
            obj.elements = mFEM.Vector(length(X));
            
            spmd 
                x = getLocalPart(X);
                id = globalIndices(X,2);
                
                elem = cell(length(x),1);
                
                % Loop through the grid, creating elements for each cell
                for i = 1:length(x);     
                    elem{i} = feval(['mFEM.elements.', elem_type], id(i), x(:,i));    
                end
                
                
                 codist = codistributor1d(codistributor1d.unsetDimension, ...
                    codistributor1d.unsetPartition, global_size);
                 elem_dist = codistributed.build(elem, codist);
            end
            obj.elements = elem_dist;
        end
    end
end

