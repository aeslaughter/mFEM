classdef FEmesh < mFEM.handle_hide
    %Class for managing and generating FEM spaces
    
    properties (SetAccess = private, GetAccess = public)
        n_elements = uint32([]);    % no. of elements in mesh
        n_dim = uint32([]);         % no. of space dimensions
        n_dof = uint32([]);         % total no. of global dofs
        n_dof_node = uint32(1);     % no. of dofs per node      
        element;                    % empty array of elements [Element]
        initialized = false;        % initialization state
        map = ...                   % struct  of node, elem, and dof maps 
            struct('node', [], 'elem', uint32([]), 'dof', uint32([]),...
            'boundary', uint32([]));       
        opt = ...                   % struct of default user options
            struct('space', 'scalar', 'type', 'CG', 'time', true);
    end
    
    properties (GetAccess = {?mFEM.System}, SetAccess = private)
        boundary_id = uint32([]);   % appended when add_boundary is called
        local_n_dim = uint32([]);   % local dimensions of elements
    end
    
    methods (Access = public)
        function obj = FEmesh(varargin)
            % FEMESH Creates a finite element mesh object.
            % 
            % Syntax:
            %   obj = FEmesh()
            %   obj = FEmesh(name)
            %   obj = FEmesh(..., 'PropertyName', PropertyValue)
            %
            % Description:
            %   obj = FEmesh() creates mesh object with default settings,
            %         equivalent to: 
            %           obj = FEmesh('Quad4','Space','scalar','Type','CG')
            %
            %   obj = FEmesh(name) allows the element type to be set
            %
            %   obj = FEmesh(..., 'PropertyName', PropertyValue) allows
            %           user to customize the behavior of the FE mesh.
            %
            % Options:
            %   'Type' = 'CG' (default) or 'DG'
            %               allows the type of element conectivity to be set
            %               as continous (CG) or discontinous (DG)
            %   'Space' = 'scalar', 'vector', <number>
            %               allows the type of FEM space to be set: scalar
            %               sets the number of dofs per node to 1, vector
            %               sets it to the no. of space dimension, and
            %               specifing a number sets it to that value.

            % Parse the user-defined options
            obj.opt = gather_user_options(obj.opt, varargin{:});
                     
            % Initialize the element property (the type does not matter)
            obj.element = mFEM.Quad4.empty;
        end
        
        function add_element(obj, type, nodes)
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
            %   >> mesh = FEmesh();
            %   >> mesh.add_element(1,'Quad4', [0,0; 1,0; 1,1; 0,1]);
            %
            %   Note, the optional space string for the element is not
            %   accepted here, as it is defined across the entire mesh
            %   when the FEmesh object is created.
            
            % Indicate that the object must be reinitialized
            obj.initialized = false;
                        
            % Create the element
            id = length(obj.element) + 1;
            obj.element(id) = feval(['mFEM.', type],...
                id, nodes, obj.opt.space);
            
            % Append the element and node maps                    
            nn = obj.element(end).n_nodes;
            obj.map.elem(end+1:end+nn,:) = id;
            obj.map.node(end+1:end+nn,:) = obj.element(id).nodes;
        end
        
        function init(obj)
            % Initialization function.
            %
            % Syntax:
            %   init()
            %
            % This function should be called before the FEmesh object is
            % used, it computes all of the necessary degree-of-freedom
            % information for the elements.
            
            % The no. of elements
            obj.n_elements = length(obj.element); 
            
            % Spatial dimension and local dimensions
            obj.n_dim = obj.element(1).n_dim;
            obj.local_n_dim = obj.element(1).local_n_dim;

            % Determine the no. of dofs per node
            if strcmpi(obj.opt.space, 'scalar');
                obj.n_dof_node = 1;
                
            elseif strcmpi(obj.opt.space, 'vector');
                obj.n_dof_node = obj.n_dim;
                
            elseif isnumeric(obj.opt.space);
                obj.n_dof_node = obj.opt.space;
                
            else
                error('FEmesh:FEmesh', 'The element space, %s, was not recongnized.',obj.opt.space);
            end
            
            % Computes the global degree-of-freedom map
            obj.compute_dof_map();
            
            % Compute the total number of dofs
            obj.n_dof = length(unique(obj.map.dof)) * obj.n_dof_node;
            
            % Locates the neighbor elements for each element
            obj.find_neighbors();

            % Set the initilization flag to true
            obj.initialized = true;
        end
                
        function add_boundary(obj, varargin)
            % Labels elements with numeric tags.
            %
            % This should be called after initialization, and id = 0
            % indicates is the default case (i.e., not identified). It is
            % also possible to specify multiple ids for an element.
            %
            % Syntax:
            %   add_boundary(loc, id)
            %   add_boundary(loc, value, id)
            %   add_boundary(id)
            %
            % Description:
            %   add_boundary(loc,id) adds a boundary id to a side;
            %   loc can be one of the following: 'left', 'right', 'top', 
            %   'bottom', 'front', or 'back'. The id is an interger that
            %   the user specifies for identifing the boundary.
            %
            %   add_boundary(loc,value,id) adds a boundary id based on
            %   the x,y,z location, so loc can be: 'x', 'y', or 'z'. The
            %   value is a location for the given coordinate that if an
            %   element contains will be specified with the given boundary
            %   id. For example, add_boundary('x',2,1) will mark all
            %   elements with a node that has a value of x = 2 with id 1.
            %
            %   add_boundary(id) all unidentified elements and sides
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
            
            % Initialize the new column of the map
            obj.map.boundary(:,end+1) = zeros(size(obj.map.elem), 'uint32');

            % Locate the elements
            idx = obj.map.node(:,col) == value;
            obj.map.boundary(idx,end) = id;
            e = unique(obj.map.elem(idx),'R2012a');
            
            % Loop through each element on the boundary
            for i = 1:length(e);

               % Add id to the element
               elem = obj.element(e(i));
               elem.on_boundary = true;
               elem.boundary_id(end+1) = id;
   
               % Locate dofs that meet critiera
               idx = elem.nodes(:,col) == value;
               dof = sort(elem.global_dof(idx));

               % Loop through sides, mark side if dofs match
               for s = 1:elem.n_sides;
                   dof_s = sort(elem.global_dof(elem.side_dof(s,:)));
                   if all(dof_s == dof)
                       elem.side(s).on_boundary = true;                       
                       elem.side(s).boundary_id = id;
                   end
               end
            end
            
            % Update mesh level boundary id storage (see get_dof for use)
            obj.boundary_id(end+1) = id;
        end

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
                        idx = idx + uint32(obj.map.boundary(:,col(c)) ~= id); 
                    
                    % Equal
                    else
                        idx = idx + uint32(obj.map.boundary(:,col(c)) == id);
                    end
                end

                % Collect the dofs
                dof = unique(obj.map.dof(logical(idx)));
                
            % Use all dofs 
            else
                dof = unique(obj.map.dof);   
            end

            % Vector FE space (uses transform_dof of an element)
            if obj.n_dof_node > 1;
                elem = obj.element(1);
                dof = elem.transform_dof(dof);
            end
        end
                
        function out = get_nodes(obj,varargin)
            % Return the spatial position of the nodes for the mesh
            %
            % Syntax:
            %   get_nodes()
            %   get_nodes(dim)
            %
            
            out = unique(obj.map.node, 'rows', 'stable'); 
            if nargin == 2 && isnumeric(varargin{1}) && varargin{1} <= obj.n_dim;
               out = out(:,varargin{1});
            end
        end
        
        function grid(obj, type, varargin)  
            % Create a mesh (1D, 2D, or 3D)
            % 
            % Syntax:
            %   grid(type, x0, x1, xn)
            %   grid(type, x0, x1, y0, y1, xn, yn)
            %   grid(type, x0, x1, y0, y1, z0, z1, xn, yn, zn)
            %
            % Input:
            %   type = name of element (e.g., 'Quad4')
            %   x0,x1 = Mesh limits in x-direction
            %   y0,y1 = Mesh limits in y-direction
            %   z0,z1 = Mesh limits in z-direction
            %   xn,yn,zn = Num. of elements in x,y,z direction
           
            % Display wait message
            if obj.opt.time;
                ticID = tmessage('Generating mesh...');
            end
            
            % Check the current element type is supported
            switch type;
                case {'Line2', 'Line3'};
                    obj.gen1Dgrid(type, varargin{:});
                case {'Quad4', 'Tri3', 'Tri6'};
                    obj.gen2Dgrid(type, varargin{:});
                otherwise
                    error('FEmesh:grid','Grid generation is not supported for the %s elem', type);
            end
            
            % Complete message
            if obj.opt.time;
                tmessage(ticID);
            end;
        end
       
        function plot(obj, varargin)
            %PLOT Plots mesh and/or data (see FEplot function for details)
            FEplot(obj,varargin{:});
        end
    end
    
    methods (Access = private)
        function compute_dof_map(obj)
            % Calculates the global degree-of-freedom map
            
            % Display wait message
            if obj.opt.time;
                ticID = tmessage('Computing the degree-of-freedom map...');
            end

            % Place the correct type of map in public property
            switch obj.opt.type;
                case 'CG'; 
                    [~,~,obj.map.dof] = unique(obj.map.node, 'rows','stable');
                case 'DG'; 
                    obj.map.dof = (1:size(obj.map.node,1))';
            end

            % Update the elements with the global dof
            for e = 1:obj.n_elements;
                elem = obj.element(e);
                elem.global_dof = obj.map.dof(obj.map.elem == e,:);
            end
            
            % Complete message
            if obj.opt.time;
                tmessage(ticID);
            end;
        end
        
        function find_neighbors(obj)
            % Locates elements that share a side
            %
            % Syntax:
            %   find_neighbors()   

            % Display wait message
            if obj.opt.time;
                ticID = tmessage('Locating neighbor elements...');
            end

            % Create a dof map for searching out neighbors
            [~, ~, CGdof] = unique(obj.map.node, 'rows','stable');
                                    
            % Loop through elements and store ids of neighboring elements
            for e = 1:obj.n_elements;
                elem = obj.element(e);      % current element

                % Define vectors of nodal values
                a = elem.nodes;

                % Compute distance of all nodes from center of this element
                cntr = mean(elem.nodes);
                d = max(obj.map.node - repmat(cntr,length(obj.map.node),1),[],2);
                
                % Create a local set of nodes to search
                L = elem.hmax();                % max distance
                cmap = obj.map.elem(d < L);     % elemnts closer than L
                cidx = cmap ~= e;               % exclude current element
                cmap = cmap(cidx);              % update cmap

                % Build list of nodes to search through 
                b = obj.map.node(d < L, :);
                b = b(cidx,:);            

                % Determine where the a and b vectors are the same
                [m,n] = size(a);
                index = zeros(length(b),n); 
                for i = 1:m;
                    aa = repmat(a(i,:), size(b,1), 1);
                    index(:,i) = all(aa == b, 2);
                end

                % Define the neighor elements
                all_neighbors = cmap(any(index,2));
                neighbors = unique(all_neighbors);

                % Locate the neighbor that shares all nodes on side
                occ = histc(all_neighbors, neighbors);
                id = neighbors(occ > obj.n_dim - 1);

                % Build global dof for sides of the current element
                dof_e = CGdof(obj.map.elem == e);
                dof_e = sort(dof_e(elem.side_dof),2);
                
                % Loop through the neighbor elements
                for n = 1:length(id);
                    
                    % Get the neighbor element
                    neighbor = obj.element(id(n));

                    % Skip if this pairing was already done
                    if any(elem.neighbors == neighbor);
                       continue;
                    end    

                    % Get the neighbor element dofs
                    dof_n = CGdof(obj.map.elem == id(n));
                    dof_n = sort(dof_n(neighbor.side_dof),2);

                    % Locate where the sides are the same
                    [~, s_e, s_n] = intersect(dof_e, dof_n, 'rows', 'R2012a');

                    % Store the values for the current element
                    elem.side(s_e).neighbor = neighbor;
                    elem.side(s_e).neighbor_side = uint32(s_n);  

                    % Store the values for the neighbor element
                    neighbor.side(s_n).neighbor = elem;
                    neighbor.side(s_n).neighbor_side = uint32(s_e);  

                    % Store current element in vector for comparision 
                    neighbor.neighbors(end+1) = elem;
                end  
            end
            
            % Complete message
            if obj.opt.time;
                tmessage(ticID);
            end;
         end 
        
        function gen1Dgrid(obj,type, x0, x1, xn)
            % Generate the 1D mesh (see grid)

            % Generate the generic grid points
            x = x0 : (x1-x0)/xn : x1;

            % Loop through the grid, creating elements for each cell
            for i = 1:length(x)-1;           
                % Define the nodes of the cell
                nodes(1) = x(i);
                nodes(2) = x(i+1);

                % Add the element(s)
                obj.add_element(type, nodes');
            end
        end
        
        function gen2Dgrid(obj, type, x0, x1, y0, y1, xn, yn)
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
                    switch type;
                        case {'Quad4'};
                            obj.add_element(type, nodes);
                            
                        case {'Tri3', 'Tri6'};
                            obj.add_element(type, nodes(1:3,:));
                            obj.add_element(type, nodes([1,3,4],:));
                    end
                end
            end
        end

        function s_n = match_side(obj, a, N)
            % Returns the side index of the neighbor that matches 
                        
            % Get the handle to the neighbor element
            if ~isa(N,'mFEM.Element');
                N = obj.element(N); % neighbor element
            end
            
            % Get the nodes for the neighbor element
            b = N.nodes;
            
            % Determine where the a and b vectors are the same
            [m,n] = size(a);
            rowIdx = 1:m;
            repmatRowIdx = rowIdx(:, ones(length(b),1));
            index = zeros(length(b),n); 
            for i = 1:m;
                 aa = a(i,:);
               % aa = repmat(a(i,:), size(b,1), 1);
                index(:,i) = all(aa(repmatRowIdx(:), 1:end) == b, 2);
            end

            % Locate the nodes that are common 
            idx = sort(find(any(index,2)));
            side = sort(N.side_dof,2);
            aa = repmat(idx',length(side),1);
            s_n = find(all(side == aa, 2));
        end
        
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
        
        function id_empty_boundary(obj, id)
            % Mark all elements sides that were not marked previously
    
            % Location for new row of boundary
            col = size(obj.map.boundary, 2) + 1;
            
            % Initialize the new row
            obj.map.boundary(:,col) = zeros(size(obj.map.elem), 'uint32');
            
            % Loop through each element on the boundary
            for e = 1:obj.n_elements;
                
                % Current element
                elem = obj.element(e);

                % Test if element is on the boundary
                if isempty(elem.on_boundary) || elem.on_boundary 
                    
                    % Loop through the sides
                    for s = 1:elem.n_sides;

                        % If side is on the boundary and contains no
                        % boundary ids, assign the desired id
                        if (elem.side(s).on_boundary || isempty(elem.side(s).neighbor)) ...
                                && isempty(elem.side(s).boundary_id);
                            
                            % Be sure to tag side as boundary
                            elem.side(s).on_boundary = true;
                            
                            % Apply the id to the element side data
                            elem.side(s).boundary_id = id;
                            
                            % Update the mesh boundary_id map
                            dof = elem.get_dof(s,'global');
                            for i = 1:length(dof);
                                idx = obj.map.dof == dof(i);
                                obj.map.boundary(idx, col) = id;
                            end
                        end
                    end
                end
            end
            
            % Update mesh level boundary matrix
            obj.boundary_id(end+1) = id;
        end         
    end
end