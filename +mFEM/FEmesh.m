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
        n_dof = [];             % total number of global degrees of freedom [double]
        initialized = false;    % initialization state [boolean]
        
        % Structure containing node, elem, and dof maps for the mesh 
        map = struct('node', [], 'elem', uint32([]), 'dof', uint32([]),...
            'boundary_id', uint32([]));
    end
    
    % Private properties
    properties (Access = private)
        CG_dof_map = uint32([]); % needed to compute neighbors for both CG and DG
        boundary_id = uint32([]); % appended each time add_boundary_id is called
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
            
            % Add the bin directory
            addpath('./bin');
            
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
            
            % Initialize the element property
            obj.element = feval(['mFEM.',obj.element_type,'.empty']);
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
            
            obj.element(id) = feval(['mFEM.', obj.element_type],...
                id, nodes, obj.space);
            
            % Append the element and node maps                    
            nn = obj.element(end).n_nodes;
            obj.map.elem(end+1:end+nn,:) = id;
            obj.map.node(end+1:end+nn,:) = obj.element(id).nodes;
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
            
            % Compute the total number of dofs
            obj.n_dof = length(unique(obj.map.dof));
            if strcmpi(obj.space,'vector');
                obj.n_dof = obj.n_dof * obj.n_dim;
            end
            
            % Locates the neighbor elements for each element
            obj.find_neighbors();
 
            % Set the initilization flag to true
            obj.initialized = true;
        end
                
        % Adds identification to elements
        function add_boundary(obj, varargin)
            % Labels elements with numeric tabs.
            %
            % This should be called after initialization, and id = 0
            % indicates is the default case (i.e., not identified). It is
            % also possible to specify multiple ids for an element.
            %
            % Syntax:
            %   add_boundary_id(loc, id)
            %   add_boundary_id(loc, value, id)
            %   add_boundary_id(id)
            %
            % Description:
            %   add_boundary_id(loc,id) adds a boundary id to a side;
            %   loc can be one of the following: 'left', 'right', 'top', 
            %   'bottom', 'front', or 'back'. The id is an interger that
            %   the user specifies for identifing the boundary.
            %
            %   add_boundary_id(loc,value,id) adds a boundary id based on
            %   the x,y,z location, so loc can be: 'x', 'y', or 'z'. The
            %   value is a location for the given coordinate that if an
            %   element contains will be specified with the given boundary
            %   id. For example, add_boundary_id('x',2,1) will mark all
            %   elements with a node that has a value of x = 2 with id 1.
            %
            %   add_boundary_id(id) all unidentified elements and sides
            %   that are on a boundary are tagged with id

            % Check if system is initialized
            if ~obj.initialized;
                error('ERROR: The mesh must be initialized');
            end
                        
            % Special case, tag all untagged
            if nargin == 2 && isnumeric(varargin{1});
                obj.id_empty_boundary(varargin{1});
                return;
            end

            % Collect the input from the user
            [col, value, id] = obj.parse_boundary_id_input(varargin{:});
            
            % Locate the elements
            idx = obj.map.node(:,col) == value;
            obj.map.boundary_id(idx,end+1) = id;
            e = unique(obj.map.elem(idx));
            
            % Loop through each element on the boundary
            for i = 1:length(e);

               % Add id to the element
               elem = obj.element(e(i));
               elem.boundary_id(end+1) = id;
   
               % Locate dofs that meet critiers
               elem_nodes = elem.nodes;
               idx = elem_nodes(:,col) == value;
               dof = sort(elem.global_dof(idx));

               % Loop through sides, mark side if dofs match
               for s = 1:elem.n_sides;
                   if all(sort(elem.side(s).global_dof,1) == dof)
                       elem.side(s).boundary_id = id;
                   end
               end
            end
            
            % Update mesh level boundary id storage (see get_dof for use)
            obj.boundary_id(end+1) = id;
        end

        % Return the global dogs for entire mesh
        function dof = get_dof(obj, varargin)
            % The global degrees of freedom, account for type of space
            %
            % Syntax:
            %   get_dof();
            %   get_dof(boundary_id)
            %   get_dof(boundary_id, 'ne');

            test = 'e';
            if nargin == 3;
                test = varargin{2};
            end

            % Limit to the specified boundary id
            if nargin >= 2 && isnumeric(varargin{1});
                % Define the id to search
                id = varargin{1};
                
                % Initialized the indices
                idx = uint32(zeros(size(obj.map.dof)));
                
                % Locate the colums of boundary_id map to search
                col = find(obj.boundary_id == id);
                
                % Loop the through the columns associated with the id and
                % update the idx
                for c = 1:length(col)
                    % Not equal
                    if strcmpi(test,'ne');
                        idx = idx + uint32(obj.map.boundary_id(:,col(c)) ~= id); 
                    
                    % Equal
                    else
                        idx = idx + uint32(obj.map.boundary_id(:,col(c)) == id);
                    end
                end

                % Collect the dofs
                dof = unique(obj.map.dof(logical(idx)));
                
            % Use all dofs 
            else
                dof = unique(obj.map.dof);   
            end

            % Vector FE space (uses transform_dof of an element)
            if strcmpi(obj.space, 'vector');
                elem = obj.element(1);
                dof = elem.transform_dof(dof);
            end
        end
        
        % Return the spatial position of nodes
        function out = get_nodes(obj,varargin)
            % Return the spatial position of the nodes for the mesh
            %
            % Syntax:
            %   get_nodes()
            %   get_nodes(dim)
            %
            
            out = unique(obj.map.node, 'rows'); 
            if nargin == 2 && isnumeric(varargin{1}) && varargin{1} <= obj.n_dim;
               out = out(:,varargin{1});
            end
        end
            
        % Tool for generating a 2D mesh
        function grid(obj, varargin)  
            % Create a mesh (1D, 2D, or 3D)
            % 
            % Syntax:
            %   grid(x0, x1, xn)
            %   grid(x0, x1, y0, y1, xn, yn)
            %   grid(x0, x1, y0, y1, z0, z1, xn, yn, zn)
            %
            % Input:
            %   x0,x1 = Mesh limits in x-direction
            %   y0,y1 = Mesh limits in y-direction
            %   z0,z1 = Mesh limits in z-direction
            %   xn,yn,zn = Num. of elements in x,y,z direction
           
            % Display wait message
            tic;
            disp('Generating mesh...'); 
            
            % Check the current element type is supported
            switch obj.element_type;
                case 'Linear2';
                    obj.gen1Dgrid(varargin{:});
                case 'Quad4';
                    obj.gen2Dgrid(varargin{:});
                otherwise
                    error('FEmesh:grid','Grid generation is not supported for the %s element', obj.element_type);
            end
            
            % Complete message
            disp(['    ...Completed in ', num2str(toc),' sec.']);
   
            % Initialize
            obj.initialize();
        end
       
        % Tool for plotting the mesh with element and node labels
        function plot(obj, varargin)
            % Plots the mesh with element and node numbers (global)
            % (currently this only works with 2D and 3D (not tested in 3D)
            %
            % Syntax:
            %   plot();
            %   plot(data);
            %   plot(...,'PropertyName',<PropertyValue>);
            %
            % Descrtiption:
            %
            % Properties:
            %   'ElementLabels' = true or false
            %   'NodeLabels' = true or false
            
            % Collect the input 
            [data, opt] = obj.parse_plot_input(varargin{:});

            % Check that mesh is initialized
            if (~obj.initialized)
                error('ERROR: FEmesh must be initialized before plot() is called.');
            end
            
            % Warning for 1D case
            if all(obj.n_dim ~= [2,3]);
               error('ERROR: FEmesh only works for 2D and 3D meshes'); 
            end
            
            % Create the figure
            if opt.newfigure;
                figure; hold on;
            end
            
            % Loop through each node
            for e = 1:obj.n_elements;
                
                % Gather nodes
                elem = obj.element(e);              % current element          
                node = num2cell(elem.nodes,1);      % col. cell of nodes
                cntr = num2cell(mean(elem.nodes,1));% mean location of nodes
                
                % Create patch
                if ~isempty(data);
                    dof = elem.get_dof();
                    patch(node{:}, data(dof));
                else
                    patch(node{:}, 'b', 'FaceColor','none');
                end
                
                % Show label
                if opt.elementlabels;
                    text(cntr{:},num2str(e),'FontSize',14,...
                        'BackgroundColor','k','FontWeight','Bold',...
                        'Color','w','HorizontalAlignment','center');
                end
                        
            end
            
            % Loop through the unique nodes
            if opt.nodelabels;
                N = unique(obj.map.node,'rows','stable');
                for i = 1:length(N);
                   node = num2cell(N(i,:));
                    text(node{:},num2str(i),'FontSize',12,'Color','w',...
                        'BackgroundColor','b','HorizontalAlignment','center');
                end
            end
        end 
    end
    
    % Private methods
    methods (Access = private)
        % Computes the dof maps
        function compute_dof_map(obj)
            % Calculates the global degree-of-freedom map
            % (this is used by both CG and DG for finding neighbors)
            
            % Display wait message
            tic;
            disp('Computing the degree-of-freedom map...');
            
            % Compute the CG dof map
            [~,~,obj.CG_dof_map] = unique(obj.map.node, 'rows','stable');

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
            
            % Complete message
            disp(['    ...Completed in ', num2str(toc),' sec.']);
            
        end
        
        % 1D mesh generation
        function gen1Dgrid(obj, x0, x1, xn)
            % Generate the 1D mesh (see grid)

            % Generate the generic grid points
            x = x0 : (x1-x0)/xn : x1;

            % Loop through the grid, creating elements for each cell
            for i = 1:length(x)-1;           
                % Define the nodes of the cell
                nodes(1) = x(i);
                nodes(2) = x(i+1);

                % Add the element(s)
                obj.add_element(nodes');
            end
        end
        
        % 2D mesh generation
        function gen2Dgrid(obj, x0, x1, y0, y1, xn, yn)
            % Generate the 2D mesh (see grid)

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
                    switch obj.element_type;
                        case {'Quad4'};
                            obj.add_element(nodes);
                            
                        case {'Tri3'};
                            obj.add_element(nodes(1:3,:));
                            obj.add_element(nodes([1,3,4],:));
                    end
                end
            end
        end
        
        % Finds element neighbors   
        function find_neighbors(obj)
            % Locates elements that share a side
            %
            % Syntax:
            %   find_neighbors()

            % Display wait message
            tic;
            disp('Locating neighbor elements...');
                        
            % Loop through elements and store ids of neighboring elements
            for e = 1:obj.n_elements;

                % Collect necessary data
                elem = obj.element(e);                  % current element
                elem_dof = elem.global_dof;             % global dof for this element
                side_dof = zeros(size(elem.side_dof));  % storage for side dof of this element
                elem.on_boundary = false;               % initilize on_boundary flag
                
                % Loop through the sides of element
                for s = 1:elem.n_sides
                    
                    % Set the local and global dof for the current side
                    dof = sort(elem.side_dof(s,:));
                    elem.side(s).dof = dof;
                    elem.side(s).global_dof = elem_dof(dof);

                    % Collects index of neighbor elements (including corners)
                    idx = zeros(size(obj.CG_dof_map));
                    for i = 1:length(side_dof(s,:));
                        idx = idx + (obj.CG_dof_map == elem.side(s).global_dof(i));
                    end
                    idx(obj.map.elem == e) = 0; %exclude current element

                    % Gathers neighbor element ids; a side can have 
                    % multiple neighbors, this is allowed for adding 
                    % adaptivity in the future
                    x = obj.map.elem(logical(idx));
                    count = zeros(size(x));
                    for i = 1:length(x);
                        count(i) = sum(x(i) == x);
                    end

                    % Identify the shared id's, corners don't count for
                    % multi-dimensional elements
                    id = unique(x(count > (obj.n_dim - 1)));
                    elem.side(s).neighbor = id;
                    
                    % Update the element data structure
                    elem.side(s).neighbor = obj.element(id);
                    
                    % Locate the matching side
                    obj.match_sides(elem,s);
                    
                    % Set the status of the neighbor flag for this element,
                    % if any side does not have a nieghbor then set the
                    % on_boundary flag to true for this elment
                    if isempty(id);
                        elem.on_boundary = true;
                        elem.side(s).on_boundary = true;
                    else
                        elem.side(s).on_boundary = false;
                    end   
                end
            end
            
            % Complete message
            disp(['    ...Completed in ', num2str(toc),' sec.']);
        end 
        
        % Match elem and neighbor side ids
        function match_sides(~, elem, s_e)
            % Returns the side index of the neighbor that matches 
            
            % Loop through the neighbors of the current side
            for i = 1:length(elem.side(s_e).neighbor);
                neighbor = elem.side(s_e).neighbor(i);
            
                % Collect the nodes
                a = elem.nodes(elem.side_dof(s_e,:));   % current element
                b = neighbor.nodes;                     % neighbor element
           
                % Locate where a & b match
                [~,~,ib] = intersect(a,b,'rows','R2012a');

                % Locate the matching side on the neighbor element
                C = b(sort(ib),:);             
                for s_n = 1:neighbor.n_sides;
                    if isequal(C, neighbor.nodes(sort(neighbor.side_dof(s_n,:)),:));
                       elem.side(s_e).neighbor_side(i) = s_n;
                       break;
                    end
                end 
            end  
        end
        
        % Sub-method for parsing boundary_id input
        function [col, value, id] = parse_boundary_id_input(obj, varargin)
            
             % The location flag
             loc = varargin{1};
             switch loc;
                 
                case {'left', 'right'};
                    if strcmpi('right', loc);
                        value = max(obj.map.node(:,1));
                    else
                        value = min(obj.map.node(:,1));
                    end
                    col = 1;
                    id = varargin{2};
                    
                case {'top', 'bottom'};
                    if obj.n_dim < 2;
                        error('ERROR: 1D case does not have a ''%s'' location.', loc);
                    end
                    if strcmpi('top', loc);
                        value = max(obj.map.node(:,2));
                    else
                        value = min(obj.map.node(:,2));
                    end
                    col = 2;
                    id = varargin{2};
                
                 case {'front', 'back'};
                    if obj.n_dim < 3;
                        error('ERROR: The ''%s'' location is only valid in 3D', loc);
                    end
                    if strcmpi('back', loc);
                        value = max(obj.map.node(:,3));
                    else
                        value = min(obj.map.node(:,3));
                    end      
                    col = 3;
                    id = varargin{1};
                
                otherwise
                    if nargin ~= 3;
                        error('ERROR: Invalid input.');
                    end
                    
                    switch lower(loc);
                        case 'x'; col = 1;
                        case 'y'; col = 2;
                        case 'z'; col = 3;
                        otherwise
                            error('ERROR: Unknown location specifier, %s', loc);
                    end
                    
                    value = varargin{2};
                    id = varargin{3};
            end  
        end
        
        % Method for tagging unmarked boundaries
        function id_empty_boundary(obj, id)
            % Mark all elements sides that were not marked previously
    
            % Location for new row of boundary_id
            col = size(obj.map.boundary_id, 2) + 1;
            
            % Loop through each element on the boundary
            for e = 1:obj.n_elements;
                
                % Current element
                elem = obj.element(e);

                % Test if element is on the boundary
                if elem.on_boundary;
                    
                    % Loop through the sides
                    for s = 1:elem.n_sides;
                        
                        % If side is on the boundary and contains no
                        % boundary ids, assign the desired id
                        if elem.side(s).on_boundary ...
                                && isempty(elem.side(s).boundary_id);
                            
                            % Apply the id to the element side data
                            elem.side(s).boundary_id = id;
                            
                            % Update the mesh boundary_id map
                            dof = elem.side(s).global_dof;
                            for i = 1:length(dof);
                                idx = obj.map.dof == dof(i);
                                obj.map.boundary_id(idx, col) = id;
                            end
                        end
                    end
                end
            end
            
            % Update mesh level boundary_id matrix
            obj.boundary_id(end+1) = id;
        end 
        
        % Function for parsing the plot input data (see plot)
        function [data, opt] = parse_plot_input(obj, varargin)
            % Function for parsing input data to plot function
            %
            % See the help for the plot function for input detals
            
            % Assume empty data  input
            data = [];
            
            % Define user properties
            opt.elementlabels = true;
            opt.nodelabels = true;
            opt.newfigure = true;

            % Parse the input
            if nargin >= 2;
                
                % The first input is always the data
                data = varargin{1};
                
                % Check that the data is sized correctly
                if ~isnumeric(data) && length(data) ~= obj.n_dof;
                    error('FEmesh:plot', 'Data not formated correctly, it must be a vector of numbers of length %d.', obj.n_dof);
                end

                % When using data, disable the labels by default and do not
                % create a new figure with each call
                opt.elementlabels = false;
                opt.nodelabels = false;
                opt.newfigure = false;
                
                % Collect options supplied by the user
                if nargin >= 2;
                    opt = gatheruseroptions(opt,varargin{2:end});
                end
            end
        end
    end
end