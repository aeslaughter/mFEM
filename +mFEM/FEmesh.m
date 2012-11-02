classdef FEmesh < mFEM.handle_hide
    %FEMESH Class for managing and generating FEM spaces.
    % This class handles all mesh and degree-of-freedom operations for
    % implementing the finite element method.

    properties (SetAccess = private, GetAccess = public)
        n_elements = uint32([]);    % no. of elements in mesh
        n_dim = uint32([]);         % no. of space dimensions
        n_dof = uint32([]);         % total no. of global dofs
        n_dof_node = uint32(1);     % no. of dofs per node      
        element;                    % empty array of elements
        initialized = false;        % initialization state
        boundary_id = uint32([]);   % list of boundary ids
        map = ...                   % structure of node, elem, dof, and boundary maps 
            struct('node', [], 'elem', uint32([]), 'dof', uint32([]),...
            'boundary', uint32([]));       
        opt = ...                   % struct of default user options
            struct('space', 'scalar', 'type', 'CG', 'time', true);
    end
    
    properties (Access = private, GetAccess = {?mFEM.System})
        local_n_dim = uint32([]);   % local dimensions of elements (see BUILD_SIDE)
    end
    
    methods (Access = public)
        function obj = FEmesh(varargin)
            % FEMESH Creates a finite element mesh object.
            % 
            % Syntax
            %   obj = FEmesh()
            %   obj = FEmesh('PropertyName', PropertyValue)
            %
            % Description
            %   obj = FEmesh() creates mesh object with default settings,
            %         equivalent to: 
            %         obj = FEmesh('Space','scalar','Type','CG')
            %
            %   obj = FEmesh('PropertyName', PropertyValue) allows
            %           user to customize the behavior of the FE mesh.
            %
            % FEMESH Property Descriptions
            %   Type
            %       {'CG'} | 'DG'
            %       Allows the type of element conectivity to be set
            %       as continous (CG) or discontinous (DG)
            %   Space
            %       {'scalar'} | 'vector' | numeric
            %       Allows the type of FEM space to be set: scalar sets the 
            %       number of dofs per node to 1, vector sets it to the no. 
            %       of space dimension, and specifing a number sets it to 
            %       that value.

            % Parse the user-defined options
            obj.opt = gather_user_options(obj.opt, varargin{:});
                     
            % Initialize the element property (the type does not matter)
            obj.element = mFEM.elements.Quad4.empty;
        end
        
        function grid(obj, type, varargin)
            %GRID Create a mesh (1D, 2D, or 3D)
            % 
            % Syntax
            %   grid(type, x0, x1, xn)
            %   grid(type, x0, x1, y0, y1, xn, yn)
            %   grid(type, x0, x1, y0, y1, z0, z1, xn, yn, zn)
            %
            % Description
            %   grid(...) creates a grid across the prescribed limits with
            %   the prescribed number of elements within the limits, see
            %   input descriptions below for additonal details.
            %
            % Description of inputs
            %   type = name of element (e.g., 'Line2')
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
            obj.element(id) = feval(['mFEM.elements.', type],...
                id, nodes, obj.opt.space);
        end

        function init(obj)
            %INIT Initialization function.
            %
            % Syntax
            %   init()
            %
            % Description
            %   init() initialize the FEMESH object instance for operation,
            %   it must be called before the FEmesh object is used (except
            %   for the GRID and ADD_ELMENENT methods), it computes all of 
            %   the necessary degree-of-freedom information for the 
            %   elements.
            
            % The no. of elements
            obj.n_elements = length(obj.element);             
            
            % Build element and node maps     
            for e = 1:obj.n_elements;
                elem = obj.element(e);
                nn = elem.n_nodes;
                obj.map.elem(end+1:end+nn,:) = elem.id;
                obj.map.node(end+1:end+nn,:) = elem.nodes;
            end
       
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
                
        function add_boundary(obj, id, varargin)
            %ADD_BOUNDARY Labels elements with numeric tags.
            %
            % Syntax
            %   add_boundary(id)
            %   add_boundary(id,Limit1,...)
            %
            % Description
            %   add_boundary(id) labels all unidentified elements and sides
            %       with the specified id, which must be an interger
            %       value.
            %
            %   add_boundary(id, Limit1,...) allows 
            %       for customization to what boundaries are tagged, see 
            %       the descriptions below. It is possible to
            %       supply multiple limits For example,
            %           add_boundary(1,'left','right','x==1')
            %       will loop through each value and apply the boundary
            %       id to each.
            %
            % Descriptions of Limits
            %   Locations:
            %       'left' | 'right' | 'top' | 'bottom' | 'front' | 'back'
            %       Adds a boundary id to a side identified by a string
            %       that describes its location. Note, in 1D use 'left' or
            %       'right', in 2D 'top' and 'bottom' are added, and in 3D
            %       all options are available.
            %
            %   Functions:
            %       string | cell array of strings
            %       It is possible to tag boundaries using test functions.
            %       For example, 'x==1' will tag boundaries that have a
            %       node location at 1. Mutliple conditions should be
            %       expressed in a cell array. For example, {'x==1','y<2'}
            %       will id elements with nodes at x = 1 AND y < 2.
            %
            % Examples
            %   The following labels the right hand side of the mesh and
            %   the node point (1,1) with id 1.
            %       add_boundary(1,'right',{'x==1','y==1'});
            %
            %   The following labels the right hand side of the mesh and
            %   the elements with nodes on x == 1 OR y == 1 to 1.
            %       add_boundary(1,'right','x==1','y==1'});
            %       

            % Check that the id is a numeric value
            if ~isnumeric(id);
                error('FEmesh:add_boundary',...
                    'The boundary id must be a number');
            end
            
            % Check if system is initialized
            if ~obj.initialized;
                error('FEmesh:add_boundary',...
                    'The FEmesh object must be initialized');
            end
                        
            % Special case, tag all untagged
            if nargin == 2;
                obj.id_empty_boundary(id);
                return;
            end
                        
            % Define array of locations to test
            cstr = {'left','right','top','bottom','front','back'};
            
            % Loop through the supplied options
            for i = 1:length(varargin);
                
                % Apply id base on location flag
                if ischar(varargin{i}) && any(strcmpi(varargin{i},cstr));
                    obj.add_boundary_location(id, varargin{i});
                
                % Apply id based on function flag
                else
                    obj.add_boundary_function(id, varargin{i});
                end
            end
        end

        function dof = get_dof(obj, varargin)
            %GET_DOF Returns the global degrees of freedom.
            %
            % Syntax
            %   dof = get_dof()
            %   dof = get_dof('PropertyName',PropertyValue);
            %
            % Description
            %   dof = get_dof() returns the global degrees of freedom for
            %   the entire finite element mesh.
            %
            %   dof = get_dof('PropertyName',PropertyValue) returns the
            %   global degrees of freedom for portions of the mesh
            %   depending on the properties (see descriptions below).
            %
            % GET_DOF Property Descriptions
            %   Component
            %       scalar | 'x' | 'y' | 'z'
            %       Returns the dof associated with the vector space
            %
            %   Boundary
            %       scalar
            %       Extract the dofs for the boundary specified, where the
            %       scalar value is the numeric id added using the
            %       ADD_BOUNDARY method.
            %
            %   Index
            %       true | {false}
            %       Toggles the type of output, by default the GET_DOF
            %       returns the dofs as a logical array equal to the length
            %       of the number of degrees of freedom. However, if the
            %       actual numeric indice is desired set this property to
            %       true. You can also use the flag style input, the
            %       following are equivlent.
            %           get_dof('index',true)
            %           get_dof('-index')
            %

            % Set the default values for the varius inputs
            options.component = [];
            options.boundary = [];
            options.index = true;
            options = gather_user_options(options, varargin{:});
            
            % Return dofs for the specified boundary id
            if ~isempty(options.boundary);
                
                % Extract the column of boundary_id map
                col = obj.boundary_id == options.boundary;
                idx = obj.map.boundary(:,col);

                % Collect the dofs
                dof = unique(obj.map.dof(logical(idx)));
                
            % Return all the dofs    
            else 
                dof = unique(obj.map.dof);  
            end

            % Vector FE space (uses transform_dof of an element)
            if obj.n_dof_node > 1;
                elem = obj.element(1);
                dof = elem.transform_dof(dof);
            end

            % Extract boundaries for certain component of vector space
            if ~isempty(options.component);
                
                % Convert string to numeric
                if ischar(options.component);
                    switch options.component;
                        case 'x'; ix = 1;
                        case 'y'; ix = 2;
                        case 'z'; ix = 3;
                    end  
                else
                    ix = options.component;
                end

                % Test the dimensionality
                if ix > obj.n_dof_node;
                    error('FEmesh:get_dof','Desired vector component does not exist.');
                end
                
                % Extract the dofs
                ix = ix:obj.n_dof_node:length(dof);
                dof = dof(ix);
            end   

            % Convert to indices, if desired
            if ~options.index
                index = false(obj.n_dof, 1);
                index(dof) = true;
                dof = index;
            end
        end
                
        function out = get_nodes(obj,varargin)
            %GET_NODES Returns spatial position of the nodes for the mesh
            %
            % Syntax
            %   x = get_nodes()
            %   x = get_nodes(dim)
            %
            % Description
            %   x = get_nodes()
            %   x = get_nodes(dim)
            
            out = unique(obj.map.node, 'rows', 'stable'); 
            if nargin == 2 && isnumeric(varargin{1}) && varargin{1} <= obj.n_dim;
               out = out(:,varargin{1});
            end
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
        
        function add_boundary_location(obj, id, loc)
           %ADD_BOUNDARY_LOCATION Adds boundary id based on location flag
           %
           % Syntax
           %    add_boundary_location(id, SideString)
           %    add_boundary_location(id, SideCell)
           %
           % Description
           %
           % See Also ADD_BOUNDARY
           
           % Convert character input into a cell
           if ischar(loc)
               loc = {loc};
           end
           
           % Loop throug each location specified and apply the id
           for i = 1:length(loc);
               
               % Locate the column in the node positions and the value to
               % search for when applying the id
               [col, value] = obj.parse_boundary_location_input(loc{i});
               
               % Build a string of the condition
               func = [col,'==',num2str(value)];
               
               % Call the function based boundary function
               obj.add_boundary_function(id, func);
           end

        end
        
        function add_boundary_function(obj, id, func)
            %ADD_BOUNDARY_FUNCTION Adds boundary id based on function
            %
            % Syntax
            %    add_boundary_function(id, FuncString)
            %    add_boundary_function(id, FuncCell)
            %
            % Description
            %
            % See Also ADD_BOUNDARY
           
            % Convert character input into a cell
            if ischar(func)
               func = {func};
            end

            % Loop throug each function specified and apply the id
            for i = 1:length(func);

                % Extract the current function
                f = strtrim(func{i});

                % Determine the column to operate on, also test that the
                % mesh is the proper dimension for the desired test
                if strcmpi(f(1),'x'); 
                   col = 1;
                elseif strcmpi(f(1),'y') && obj.n_dim > 1;
                   col = 2;
                elseif strcmpi(f(1),'z') && obj.n_dim == 3;
                   col = 3;
                else
                    error('FEmesh:add_boundary_function',...
                        'The input boundary function %s is invalid.', func{i});
                end

                % Build the strings to evaluate
                str = ['obj.map.node(:,', num2str(col), ')', f(2:end)];
                str_elem = ['elem.nodes(:,', num2str(col), ')', f(2:end)];

                % Locate the nodes on the boundary
                idx(:,i) = eval(str);
            end
           
            % Combine the indices
            idx = all(idx,2);
           
            % Meld boundary map columns together if id already exists
            c = find(obj.boundary_id == id);
            if ~isempty(c);
                idx_old = obj.map.boundary(:,c);
                obj.map.boundary(:,c) = any([idx_old,idx],2);

            % Create a new column in the boundary map if id is new    
            else
                obj.map.boundary(:, end+1) = idx;
                obj.boundary_id(end+1) = id;
            end
           
            % Locate the elements on the boundary
            E = unique(obj.map.elem(idx),'R2012a');

            % Loop through each element on the boundary
            for e = 1:length(E);

               % Add id to the element
               elem = obj.element(E(e));
               elem.on_boundary = true;
               elem.boundary_id(end+1) = id;

               % Locate dofs that meet critiera
               e_idx = eval(str_elem);
               dof = sort(elem.global_dof(e_idx));

               % Loop through sides, mark side if dofs match
               for s = 1:elem.n_sides;
                   dof_s = sort(elem.global_dof(elem.side_dof(s,:)));
                   if all(dof_s == dof)
                       elem.side(s).on_boundary = true;                       
                       elem.side(s).boundary_id = id;
                   end
               end
            end
        end
        
        function [col, value] = parse_boundary_location_input(obj, loc)
            %PARSE_BOUNDARY_LOCATION_INPUT
            %
            % Syntax
            %
            % Description
            %
            % 
            
            switch loc;
                case {'left', 'right'};
                    if strcmpi('right', loc);
                        value = max(obj.map.node(:,1));
                    else
                        value = min(obj.map.node(:,1));
                    end
                    col = 'x';

                case {'top', 'bottom'};
                    if obj.n_dim < 2;
                        error('ERROR: 1D case does not have a ''%s'' location.', loc);
                    end
                    if strcmpi('top', loc);
                        value = max(obj.map.node(:,2));
                    else
                        value = min(obj.map.node(:,2));
                    end
                    col = 'y';

                 case {'front', 'back'};
                    if obj.n_dim < 3;
                        error('ERROR: The ''%s'' location is only valid in 3D', loc);
                    end
                    if strcmpi('back', loc);
                        value = max(obj.map.node(:,3));
                    else
                        value = min(obj.map.node(:,3));
                    end      
                    col = 'z';

                 otherwise
                    error('FEmesh:parse_boundary_location_input',...
                    'Unknown location specifier, %s', loc);
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