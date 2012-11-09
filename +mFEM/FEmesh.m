classdef FEmesh < mFEM.handle_hide
    %FEMESH Class for managing and generating FEM spaces.
    % This class handles all mesh and degree-of-freedom operations for
    % implementing the finite element method.
    %
    %----------------------------------------------------------------------
    % Copyright 2012 Andrew E. Slaughter
    % This software is for educational purposes only and may not be used
    % without written permession.
    %----------------------------------------------------------------------

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
            struct('space', 'scalar', 'type', 'CG', 'time', true,...
                'mixed', false);
    end
    
    properties (SetAccess = private, GetAccess = ?mFEM.System)
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
            %
            %   Space
            %       {'scalar'} | 'vector' | numeric
            %       Allows the type of FEM space to be set: scalar sets the 
            %       number of dofs per node to 1, vector sets it to the no. 
            %       of space dimension, and specifing a number sets it to 
            %       that value.
            %
            %   Mixed
            %       true | {false}
            %       Set this flag to true if the mesh contains mixed
            %       elements, this changes how the quadrature rules are
            %       applied. If true the quadrature rules are called from
            %       every element when assembled using the System class,
            %       which takes more time. The default uses the quadrature
            %       rule from the first element for all the elements.

            % Parse the user-defined options
            obj.opt = gather_user_options(obj.opt, varargin{:});
                     
            % Initialize the element property (the type does not matter)
            obj.element = mFEM.elements.Line2.empty;
        end
        
        function grid(obj, type, varargin)
            %GRID Create a mesh (1D, 2D, or 3D)
            % 
            % Syntax
            %   grid(type, x0, x1, xn)
            %   grid(type, x0, x1, y0, y1, xn, yn)
            %   grid(type, x0, x1, y0, y1, z0, z1, xn, yn, zn)
            %   grid(..., 'PropertyName', PropertyValue,...)
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
            %
            % GRID Property Descriptions
            %   pol2cart
            %       true | {false}
            %       Converts the inputted polar cordinates to cartesian
            %       coordinates when creating the elements. Only available
            %       for 2D and 3D grids. The flag style input may also be 
            %       used (i.e., '-pol2car').
           
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
            
            % Tag all boundaries as 0
            obj.id_empty_boundary(0);          
            
            % Locates the neighbor elements for each element
            obj.find_neighbors();

            % Set the initilization flag to true
            obj.initialized = true;
        end
                
        function add_boundary(obj, id, varargin)
            %ADD_BOUNDARY Labels elements with numeric tags.
            %
            % Syntax
            %   add_boundary(id,Limit1,...)
            %   add_boundary(id)
            %
            % Description
            %   add_boundary(id) labels all unidentified elements and sides
            %       with the specified id, which must be an interger
            %       value greater than zero.
            
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
            if ~isnumeric(id) || id == 0;
                error('FEmesh:add_boundary',...
                    'The boundary id must be a number greater than 0');
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
            options.index = false;
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
                dof = transform_dof(dof,obj.n_dof_node);
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

            % Convert to indices (the default), if desired 
            if ~options.index
                index = false(obj.n_dof, 1);
                index(dof) = true;
                dof = index;
            end
        end
                
        function varargout = get_nodes(obj,varargin)
            %GET_NODES Returns spatial position of the nodes for the mesh
            %
            % Syntax
            %   x = get_nodes()
            %   [x,y,z] = get_nodes()
            %   x = get_nodes(dim)
            %
            % Description
            %   x = get_nodes() returns the nodal coordinates for the
            %   entire mesh.
            %
            %   [x,y,z] = get_nodes() returns the nodal coordinates for
            %   each spatial dimension as a seperate output.
            %
            %   x = get_nodes(dim) returns the nodeal coordinates for the
            %   specified dimension, it may be a string ('x','y', or 'z')
            %   or a numeric value (1,2, or 3).
            
            % Extract all the nodes
            out = unique(obj.map.node, 'rows', 'stable');       

            % Extract only a portion of the nodes
            if nargin == 2;
                % Determine the column, if string is supplied
                if ischar(varargin{1});
                    switch lower(varargin{1});
                        case 'x'; c = 1;
                        case 'y'; c = 2;
                        case 'z'; c = 2;
                        otherwise
                            error('FEmesh:get_nodes','Input error, the string %s was not understood. Must specify ''x'', ''y'', or ''z''');
                    end
                  
                % Extract a numeric entry    
                elseif isnumeric(varargin{1});
                    c = varargin{1};
                else
                    error('FEmesh:get_nodes','Input error.');    
                end   
                
                % Test that the coordinate dimension exists
                if c > obj.n_dim;
                    error('FEmesh:get_nodes','Input error, the desired nodal coordinate does not exist.');
                end

                % Extract the component
                out = out(:,c);
            end
            
            % Determine output format
            if nargout == 1;
                varargout{1} = out;
            else
                varargout = num2cell(out,1);
            end
            
        end

        function plot(obj, varargin)
            %PLOT Plots mesh and/or data 
            % This is a wrapper on an independant function, FEPLOT,
            % that is located in the bin directory. For complete details
            % see FEPLOT.
            %
            % Syntax
            %   plot()
            %   plot(data)
            %   plot(...,'PropertyName',PropertyValue,...)
            %
            % Description
            %   plot() creates a plot of the mesh and includes that
            %   includes the element and degree of freedom labels
            %
            %   plot(data) creates a plot of the data supplied, the data
            %   must be a column vector of with the number of values equal
            %   to the number of degrees of freedom for the mesh.
            %
            %   plot(...,'PropertyName',PropertyValue,...) allows for the
            %   plot created to be customised, for a list of properties see
            %   the help for the FEPLOT function.
            %       doc FEPlot;
            %
            % See Also FEPLOT
            FEplot(obj,varargin{:});
        end
    end
    
    methods (Hidden = true, Access = private)        
        function compute_dof_map(obj)
            %COMPUTE_DOF_MAP Calculates the global degree-of-freedom map
            %
            % Syntax
            %   compute_dof_map()
            %
            % Description
            %   compute_dof_map() builds a degree of freedom map for the
            %   finite element mesh.
            
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
            % Syntax
            %   find_neighbors()   
            %
            % Description
            %   find_neighbors() locates neighboring elements and
            %   neighboring sides for the finite element mesh, it sets the
            %   various neighbor and side related properties for the
            %   elements in the mesh.

            % Display wait message
            if obj.opt.time;
                ticID = tmessage('Locating neighbor elements...');
            end

            % Create a dof map for searching out neighbors
            [~, ~, CGdof] = unique(obj.map.node, 'rows','stable');
            
            % Loop through elements and store ids of neighboring elements
            parfor e = 1:obj.n_elements;
                
                % The current element
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
                    elem.side(s_e).on_boundary = false;

                    % Store the values for the neighbor element
                    neighbor.side(s_n).neighbor = elem;
                    neighbor.side(s_n).neighbor_side = uint32(s_e);  
                    neighbor.side(s_n).on_boundary = false;

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
            %GEN1DGRID Generate the 1D mesh (see grid)
            %
            % Syntax
            %   gen1Dgrid(type, x0, x1, xn)
            %
            % Description
            %   gen1Dgrid(type, x0, x1, xn) creates a 1D grid ranging from
            %   x0 to x1 with xn number of elements.

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
        
        function gen2Dgrid(obj, type, x0, x1, y0, y1, xn, yn, varargin)
            %GEN2DGRID Generate the 2D mesh (see grid)
            %
            % Syntax
            %   gen2Dgrid(type, x0, x1, y0, y1, xn, yn)
            %
            % Description
            %   gen2Dgrid(type, x0, x1, y0, y1, xn, yn) creates a 2D grid 
            %   ranging from x0 to x1 with xn number of elements in the
            %   x-direction, similarily in the y direction.
            %
            % GRID Property Descriptions
            %   pol2cart
            %       true | {false}
            %       Converts the inputted polar cordinates to cartesian
            %       coordinates when creating the elements. The flag style 
            %       input may also be used (i.e., '-pol2car').
            
            % Collect user options
            options.pol2cart = false;
            options = gather_user_options(options,varargin{:});
            
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
                    
                    % Adust for polar input
                    if options.pol2cart;
                        [nodes(:,1),nodes(:,2)] = ...
                            pol2cart(nodes([1,4,3,2],1), nodes([1,4,3,2],2));
                    end
                    
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
        
        function add_boundary_location(obj, id, loc)
           %ADD_BOUNDARY_LOCATION Adds boundary id based on location flag
           %
           % Syntax
           %    add_boundary_location(id, SideString)
           %
           % Description
           %    add_boundary_location(id, SideString) adds a boundary id
           %    specified by the SideString value, see the help for
           %    ADD_BOUNDARY for a list of the available strings.
            %
           % See Also ADD_BOUNDARY
           
           % Locate the column in the node positions and the value to
           % search for when applying the id
           [col, value] = obj.parse_boundary_location_input(loc);

           % Build a string of the condition
           func = [col,'==',num2str(value)];

           % Call the function based boundary function
           obj.add_boundary_function(id, func);
        end
        
        function add_boundary_function(obj, id, func)
            %ADD_BOUNDARY_FUNCTION Adds boundary id based on function
            %
            % Syntax
            %    add_boundary_function(id, FuncString)
            %    add_boundary_function(id, FuncCell)
            %
            % Description
            %    add_boundary_function(id, FuncString) adds a boundary id
            %    based on the string expression in FuncString, see
            %    ADD_BOUNDARY for details.
            %
            %    add_boundary_function(id, FuncCell) adds a boundary id
            %    based on all of the string expressions in FuncCell, see
            %    ADD_BOUNDARY for details.
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
           
            % Remove the new boundary from 0 id
            obj.map.boundary(idx,1) = false;
            
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
            %   [col, value] = parse_boundary_location_input(loc)
            %
            % Description
            %   [col, value] = parse_boundary_location_input(loc) returns
            %   the column and value for boundary location strings, see
            %   ADD_BOUNDARY for details.
            %
            % See Also ADD_BOUNDARY ADD_BOUNDARY_LOCATION

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
            %ID_EMPTY_BOUNDARY Mark all unmarked boundaries
            %
            % Syntax
            %   id_empty_boundary(id)
            %
            % Description
            %   id_empty_boundary(id) tags all boundaries that are unmarked
            %   with the specified id, see ADD_BOUNDARY for details.
            %
            % See Also ADD_BOUNDARY
    
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
                                && (isempty(elem.side(s).boundary_id) || elem.side(s).boundary_id == 0);
                            
                            % Be sure to tag side and elemen as boundary
                            elem.side(s).on_boundary = true;
                            elem.on_boundary = true;
                            
                            % Apply the id to the element side data
                            elem.side(s).boundary_id = id;
                            
                            % Update the mesh boundary_id map
                            dof = elem.get_dof('Side', s); % global
                            for i = 1:length(dof);
                                idx = obj.map.dof == dof(i);  
                                obj.map.boundary(idx, col) = true;
                                
                                % Remove from unmarked boundaries
                                if col > 1
                                    obj.map.boundary(idx,1) = false;
                                end
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