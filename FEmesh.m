classdef FEmesh < handle
    % Class for managing and generating FEM spaces
    
    % Read only properties
    properties (SetAccess = private, GetAccess = public)
        n_elements = [];        % no. of elements in mesh [double]
        element;                % empty array of elements [Element]
        element_type = 'Quad4'; % element type (default is Quad4) [string]
        n_dim = [];             % no. of spatial dimensions [double]
        type = 'CG';            % FE type ('CG' or 'DG') [string]
        space = 'scalar';       % 'scalar' or 'vector' space [string]
        initialized = false;    % initialization state [boolean]
        
        % Structure containing node, elem, and dof maps for the mesh 
        map = struct('node', [], 'elem', uint32([]), 'dof', uint32([]));
    end
    
    % Private properties
    properties (Access = private)
        CG_dof_map = uint32([]); % needed to compute neighbors for both CG and DG
    end
    
    % Public methods
    methods (Access = public)
        % Constructor
        function obj = FEmesh(varargin)
            % Class constructor
            % 
            % Syntax:
            %   obj = FEmesh()
            %   obj = FEmesh(name)
            %   obj = FEmesh(name, space)
            %   obj = FEmesh(name, space, type)
            %
            % Description:
            %   obj = FEmesh() creates mesh object with default settings,
            %         equivalent: obj = FEmesh('Quad4','scalar','CG')
            %   obj = FEmesh(name) allows the element type to be set
            %   obj = FEmesh(name, space) allows the type of FEM space
            %         to be set: 'scalar' or 'vector'
            %   obj = FEmesh(name, space, type) allows the type of element
            %         conectivity to be set: 'CG' or 'DG'
            
            % Assign element name
            if nargin >= 1;
                obj.element_type = varargin{1};
            end

            % User provided the FE space ('scalar'(default) or 'vector')
            if nargin >= 2;
                obj.space = varargin{2};
            end
            
            % User provided the FE type ('CG'(default) of 'DG')
            if nargin == 3;
                obj.type = varargin{3};
            end
            
            
            obj.element = feval([obj.element_type,'.empty']);
        end
        
        % Adds an element to the mesh
        function add_element(obj, nodes)
            % Add element to FEmesh object (must be reinitialized)
            %
            % Syntax:
            %   add_element(nodes);
            %
            % Input:
            %   The input should take the exact form as the input for the
            %   element that was established when the FEmesh object was
            %   created. For example:
            %
            %   >> mesh = FEmesh('Quad4');
            %   >> mesh.add_element(1,[0,0; 1,0; 1,1; 0,1]);
            %
            %   Note, the optional space string for the element is not
            %   accepted here, as it is defined across the entire mesh
            %   when the FEmesh object is created.
            
            % Indicate that the object must be reinitialized
            obj.initialized = false;
            
            % Add the element(s)
            id = length(obj.element) + 1;
            obj.element(id) = feval(obj.element_type, id, nodes, obj.space);
            
            % Append the element and node maps                    
            nn = obj.element(end).n_nodes;
            obj.map.elem(end+1:end+nn,:) = id;
            obj.map.node(end+1:end+nn,:) = nodes;
        end
        
        % Initializes by computing degree-of-freedom maps
        function initialize(obj)
            % Initialization function.
            %
            % Syntax:
            %   initialize()
            %
            % This function should be called before the FEmesh object is
            % used, it computes all of the necessary degree-of-freedom
            % information for the elements.
            
            % The no. of elements
            obj.n_elements = length(obj.element); 
            
            % Spatial dimension
            obj.n_dim = obj.element(1).n_dim;
            
            % Computes the global degree-of-freedom map
            obj.compute_dof_map();
            
            % Locates the neighbor elements for each element
            obj.find_neighbors();
            
            % Set the initilization flag to true
            obj.initialized = true;
        end
                
        % Tool for generating a 2D mesh
        function gen2D(obj, x0, x1, y0, y1, xn, yn)  
            % Create a 2D mesh
            % 
            % Syntax:
            %   gen2d(x0, x1, y0, y1, xn, yn)
            %
            % Input:
            %   x0,x1 = Mesh limits in x-direction
            %   y0,y1 = Mesh limits in y-direction
            %   xn,yn = Num. of elements in x,y direction
            
            % Check the current element type is supported
            switch obj.element_type;
                case 'Quad4';
                    obj.gen2Dmesh(x0, x1, y0, y1, xn, yn);
                otherwise
                    error('Mesh.gen2D is not supported for this element');
            end
        end
       
        % Tool for plotting the mesh with element and node labels
        function plot(obj)
            % Plots the mesh with element and node numbers (global)
            % (currently this only works with 2D and 3D (not tested in 3D)
            
            % Warning for 1D case
            if all(obj.n_dim ~= [2,3]);
               error('ERROR: FEmesh only works for 2D and 3D meshes'); 
            end
            
            % Create the figure
            figure; hold on;
            
            % Loop through each node
            for e = 1:obj.n_elements;
                
                % Gather nodes
                elem = obj.element(e);              % current element          
                node = num2cell(elem.nodes,1);      % col. cell of nodes
                cntr = num2cell(mean(elem.nodes,1));% mean location of nodes
                
                % Create patch
                patch(node{:}, 'b', 'FaceColor','none');
                
                % Show label
                text(cntr{:},num2str(e),'FontSize',14,'BackgroundColor','k',...
                    'FontWeight','Bold','Color','w',...
                    'HorizontalAlignment','center');
            end
            
            % Loop through the unique nodes
            N = unique(obj.map.node,'rows');
            for i = 1:length(N);
               node = num2cell(N(i,:));
                text(node{:},num2str(i),'FontSize',12,'Color','w',...
                    'BackgroundColor','b','HorizontalAlignment','center');
                 
            end
        end 
    end
    
    % Private methods
    methods (Access = private)
        % Computes the dof maps
        function compute_dof_map(obj)
            % Calculates the global degree-of-freedom map
            % (this is used by both CG and DG for finding neighbors)
            
            % Compute the CG dof map
            [~,~,obj.CG_dof_map] = unique(obj.map.node, 'rows');
            
            % Place the correct type of map in public property
            switch obj.type;
                case 'CG'; 
                    obj.map.dof = obj.CG_dof_map;
                case 'DG'; 
                    obj.map.dof = (1:size(obj.map.node,1))';
            end

            % Update the elements with the global dof
            for e = 1:obj.n_elements;
                elem = obj.element(e);
                elem.global_dof = obj.CG_dof_map(obj.map.elem == e,:);
            end
        end
        
        % 2d mesh generation
        function gen2Dmesh(obj, x0, x1, y0, y1, xn, yn)
            % Generate the 2D mesh (see gen2D)

            % Generate the generic grid points
            x = x0 : (x1-x0)/xn : x1;
            y = y0 : (y1-y0)/yn : y1;

            % Loop through the grid, creating elements for each cell
            for i = 1:length(x)-1;           
                for j = 1:length(y)-1;
                    % Define the nodes of the cell
                    nodes(1,:) = [x(i), y(j)];
                    nodes(2,:) = [x(i+1), y(j)];
                    nodes(3,:) = [x(i+1), y(j+1)];
                    nodes(4,:) = [x(i), y(j+1)];
                   
                    % Add the element(s)
                    obj.add_element(nodes);
                end
            end
            
            % Initialize
            obj.initialize();
        end 
                     
        % Finds element neighbors
        function find_neighbors(obj)
            % Locates elements that share a side
           
            % Loop through elements and store ids of neighboring elements
            for e = 1:obj.n_elements;
                elem = obj.element(e); % current element
                elem_dof  = elem.global_dof; % global dof for this element
                side_dof = zeros(size(elem.side_dof)); % storage for side dof of this element
                
                % Loop through the sides of element
                for s = 1:elem.n_sides
                    
                    % Append to the global dof for the current side
                    side_dof(s,:) = elem_dof(elem.side_dof(s,:));

                    % Collects index of neighbor elements (including corners)
                    idx = zeros(size(obj.CG_dof_map));
                    for i = 1:length(side_dof(s,:));
                        idx = idx + (obj.CG_dof_map == side_dof(s,i));
                    end
                    idx(obj.map.elem == e) = 0; %exclude current element

                    % Gathers neighbor element ids, only include as a neighbor
                    % if it shares multiple nodes (i.e., corners don't count)
                    % Also, a side can have multiple neighbors, this is
                    % allowed for adding adaptivity in the future
                    x = obj.map.elem(logical(idx));
                    count = zeros(size(x));
                    for i = 1:length(x);
                        count(i) = sum(x(i) == x);
                    end
                    elem.neighbor(1).element{s} = unique(x(count > 1));
                end
                
                 % Assign global dof for the sides
                elem.global_side_dof = side_dof;             
            end
            
           % Loop through all elements and identify which sides are shared
            for e = 1:obj.n_elements;
                elem = obj.element(e); % current element

                % Global dof for the sides of the current element, it must
                % be sorted because the node order for each side differs
                % for the element vs. its neighbors for the same dofs
                a = sort(elem.global_side_dof, 2);

                % Loop through all neighboring elements
                ne = elem.neighbor.element; 
                for s = 1:length(ne); % loop through sides
                    for i = 1:length(ne{s}); % loop through neighbor elements on current side
                        elem_n = obj.element(ne{s}(i)); % current neighbor element
                        b = sort(elem_n.global_side_dof, 2); % neighbor element global side dofs
                        [~,~,ib] = intersect(a,b,'rows','R2012a'); % fine where a and b match
                        elem.neighbor.side_index{s}(i) = ib; % update current element
                    end
                end
                    

            end
        end  
    end
end