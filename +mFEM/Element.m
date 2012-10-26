classdef Element < mFEM.handle_hide & matlab.mixin.Heterogeneous
    %ELEMENT Base class for defining elements.
    % Inludes the general behavior of an element, including the node 
    % locations, id, shape functions, etc...
    %
    % This is an abstract class, as such it must be inherited to function.
    % The abstract properties and methods must be redifined in the
    % subclass, see Quad4.m for an example.
    %
    % see Quad4
    
    % Abstract Properties (must be redefined in subclass)
    properties (Abstract = true, SetAccess = protected, GetAccess = public) 
      n_sides;      % no. of sides
      side_dof;     % array of local node ids for each side
      side_type;    % defines the type of element that defines the sides
      quad;         % Instance of Gauss quadrature class to use
    end
    
    % Abstract Methods (protected)
    % (the user must redfine these in subclasse, e.g. Quad4)
    methods (Abstract, Access = protected)
        N = basis(obj, varargin)            % basis functions
        B = grad_basis(obj, varargin)       % basis function derivatives (dN/dx, ...)
        G = local_grad_basis(obj, varargin) % basis function derivatives (dN/dxi, ...)
        J = jacobian(obj, varargin)         % the Jacobian matrix for the element
    end
        
    % Public properties (read only)
    properties (SetAccess = protected, GetAccess = public)
        id = uint32([]);            % element id [double]
        n_nodes = uint32([]);       % no. of nodes [double]
        n_dim = uint32([]);         % no. of spatial dimensions
        n_dof = uint32([]);         % no. of global degrees of freedom
        n_dof_node = uint32(1);     % no. of dofs per node (scalar = 1)  
        nodes = [];                 % global coordinates (no. nodes, no. dim)
        node_plot_order;            % node plotting order (see Tri6; only needed if nodes are out of order)
    end
    
    % Public properties (read only; except FEmesh and Element)
    properties (SetAccess = {?mFEM.FEmesh, ?mFEM.Element}, SetAccess = protected, GetAccess = public)
        on_boundary;                % flag if element is on a boundary
        boundary_id = uint32([]);   % list of all boundary ids for element
        side;                       % side info, see constructor        
        local_n_dim;                % local dimensions (default is n_dim; see Truss2 for exception)
    end
    
    % Protected properties
    properties (Access = {?mFEM.FEmesh, ?mFEM.Element}, Access = protected)
       global_dof = []; % global dof for nodes of element 
       neighbors;       % storage of nieghbor elements (see FEmesh.find_neighbors)
    end    
    
    % Public Methods
    % (These methods are accessible by the user to create the element and
    % access the shape functions and other necessary parameters)
    methods (Access = public)
        function obj = Element(id, nodes, varargin)
            % Class constructor.
            %
            % This is an abstract class, it must be inherited by a subclass
            % to operate, see Quad4.m for example.
            %
            % Syntax:
            %   Element(id, nodes)
            %   Element(id, nodes, n_dof_node)
            %
            % Description:
            %   Element(id, nodes) creates an element given:
            %       id: unique identification number for this element
            %       nodes: matrix of node coordinates (global), should be 
            %              arranged as column matrix [no. nodes x no. dims]
            %
            %   Element(id, nodes, space) allows the user to
            %       customize the behavior of the element, the available
            %       properties are listed below.
            %
            % Properties:
            %   space = 'scalar', 'vector', <number>
            %               allows the type of FEM space to be set: scalar
            %               sets the number of dofs per node to 1, vector
            %               sets it to the no. of space dimension, and
            %               specifing a number sets it to that value.

            % Insert required values into object properties
            obj.id = id;
            obj.nodes = nodes;
            [obj.n_nodes, obj.n_dim] = size(nodes);
            obj.local_n_dim = obj.n_dim;

            % Change dofs per node
            if nargin == 3 && strcmpi(varargin{1},'vector');
                obj.n_dof_node = obj.n_dim;  
            elseif nargin == 3 && isnumeric(varargin{1});
            	obj.n_dof_node = varargin{1};
            end   
            
            % Determine the total number of global dofs
            obj.n_dof = obj.n_nodes * obj.n_dof_node;
            
            % Intialize the side data structure
            obj.side = struct('on_boundary', false, ...
                'boundary_id', cell(obj.n_sides,1),...
                'neighbor',[], 'neighbor_side', []);
            
            % Initialize neighbor array
            obj.neighbors = feval([class(obj),'.empty']);
        end
        
        function N = shape(obj, varargin)
            % Returns the shape functions

            % Scalar field basis functions
            N = obj.basis(varargin{:});

            % Non-scalar fields
            if obj.n_dof_node > 1;
                n = N;                          % re-assign scalar basis
                r = obj.n_dof_node;             % no. of rows
                c = obj.n_dof_node*obj.n_nodes; % no. of cols
                N = zeros(r,c);                 % size the vector basis
    
                % Loop through the rows and assign scalar basis
                for i = 1:r;
                    N(i,i:r:c) = n;
                end
            end            
        end
        
        function B = shape_deriv(obj, varargin)
            % Returns the shape function derivatives in x,y system

            % Scalar field basis functin derivatives
            B = obj.grad_basis(varargin{:});
                        
            % Non-scalar fields
            if obj.n_dof_node > 1;
                b = B;                      % Re-assign scalar basis
                r = obj.n_dof_node;         % no. of rows
                c = r*size(b,2);            % no. of cols
                B = zeros(r+1,c);           % size the vector basis

                % Loop through the rows and assign scalar basis
                for i = 1:r;
                    B(i,i:r:c)  = b(i,:);
                    B(r+1, i:r:c) = b((r+1)-i,:);
                end
            end
        end
        
        function L = hmax(obj)
            % Return maximum length between ANY two points
            L = max(obj.distance());
        end
        
        function L = hmin(obj)
            % Return maximum length between ANY two points
            L = min(obj.distance());
        end
            
        function J = detJ(obj, varargin)
            % Returns the determinate of the jacobian matrix
            J = det(obj.jacobian(varargin{:}));
        end 
        
        function varargout = get_position(obj, varargin)
            % Returns the real coordinates given xi, eta, ...
            %
            % Syntax:
            %   x = get_position(xi)
            %   [x,y] = get_position(xi,eta)
            %   [x,y,z] = get_position(xi,eta,zeta)
            %
            % Note: This function accounts for being on a side, if the
            % element is a side element then it will use the side_nodes
            % variable to give you the position of the point in the real
            % x,y,z coordinate system.
           
            % The number of spatial dimensions available
            n = obj.n_dim;
            
            % Initialize the output
            varargout = cell(n);
            
            % Loop through the dimensions and return the desired position
            for i = 1:n;
               varargout{i} = obj.shape(varargin{:})*obj.nodes(:,i); 
            end

        end
        
        function n = get_normal(obj, varargin)
            % Returns the normal vector at a given xi, eta, ...
            
            % The number of spatial dimensions available
            n = obj.n_dim; 
            
            % 2D: Side is defined by a line
            if n == 2;
                % Compute the tangent at the point
                n = obj.local_grad_basis(varargin{:}) * obj.side_nodes;
                
                % Re-arrange tangent to give the normal (outward from
                % element face is positive)
                n = [n(2), -n(1)]/norm(n);

            % 3D: Side is defined by a plane    
            elseif n == 3;
                error('Element:get_normal', 'Not yet supported');
                
            % Only defined for 1D and 2D sides
            else
                error('Element:get_normal', 'Not defined for an element with %d spatial coordinates', n);
            end 
        end
                
        function side = build_side(obj, id)
            % Build an element for the side
            
            if obj.n_dim == 3;
                error('Element:build_side','Feature not yet supported in 3D');
            end
            
            % Extract the nodes for the side
            dof = obj.side_dof(id,:);
            node = obj.nodes(dof,:);

            % Create the side element
            side = feval(['mFEM.',obj.side_type], NaN, node, obj.n_dof_node);
            
            % Set the global dofs
            side.global_dof = obj.global_dof(dof);
        end
               
        function dof = get_dof(obj, varargin)
            % The global degrees of freedom, account for type of space
            %
            % Syntax:
            %   get_dof()
            %   get_dof(s)
            %   get_dof(s,'global')
            % 
            % Description:            
            %   get_dof() returns GLOBAL dofs for the element
            %   get_dof(s) returns LOCAL dofs for the side s
            %   get_dof(s,'global') returns GLOBAL dofs for the side s
            
            % Determine the type of dof to return
            if nargin == 1; % global for entire element
                dof = obj.global_dof;
                
            elseif nargin == 2; % local for  side
                s = varargin{1};
                dof = obj.side_dof(s,:);
                
            elseif nargin == 3 && strcmpi(varargin{2},'global');
                s = varargin{1};
                dof = obj.side_dof(s,:);
                dof = obj.global_dof(dof);
            end

            % Non-scalar FE space
            if obj.n_dof_node > 1;
                dof = obj.transform_dof(dof);
            end
        end  
                        
        function D = transform_dof(obj, d)
            % Converts the dofs for vector element space
            
            n = obj.n_dof_node;         % no. of dofs per node
            D = zeros(n*length(d),1);   % size of vector space dofs
            
            % Loop through dimensions and build vector
            for i = 1:n;
                D(i:n:end) = d*n - (n-i);
            end 
        end
    end
    
    methods (Access = private)
        function d = distance(obj)
            % Compute distances between all points
           
            d = zeros(factorial(obj.n_nodes-1),1);
            k = 0;
            for i = 1:obj.n_nodes;
                for j = i+1:obj.n_nodes;
                    k = k + 1;
                    d(k) = norm(obj.nodes(i,:) - obj.nodes(j,:));
                end
            end

        end
    end
end
    