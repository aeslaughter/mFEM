classdef Element < mFEM.handle_hide & matlab.mixin.Heterogeneous
    %ELEMENT Base class for defining elements.
    % Inludes the general behavior of an element, including the node 
    % locations, id, shape functions, etc...
    %
    % This is an abstract class, as such it must be inherited to function.
    % The abstract properties and methods must be redifined in the
    % subclass, see Line2.m for an example. In general, if you need help 
    % for an element see the the help for the subclass itself.
    %
    % See Also mFEM.elements.Line2
    %
    %----------------------------------------------------------------------
    % Copyright 2012 Andrew E. Slaughter
    % This software is for educational purposes only and may not be used
    % without written permession.
    %----------------------------------------------------------------------

    % Abstract Properties (must be redefined in subclass)
    properties (Abstract = true, SetAccess = protected, GetAccess = public) 
      n_sides;      % no. of sides
      side_dof;     % array of local node ids for each side
      side_type;    % defines the type of element that defines the sides
      quad;         % Instance of Gauss quadrature class to use
    end
    
    % Abstract Methods (protected)
    % (the user must redfine these in subclasse, e.g. Line2)
    methods (Abstract, Access = protected)
        N = basis(obj, varargin)            % basis functions
        B = grad_basis(obj, varargin)       % basis function derivatives (dN/dx, ...)
        G = local_grad_basis(obj, varargin) % basis function derivatives (dN/dxi, ...)
        J = jacobian(obj, varargin)         % the Jacobian matrix for the element
    end
        
    % Public properties (read only)
    properties (SetAccess = protected, GetAccess = public)
        id = uint32([]);            % element id
        n_nodes = uint32([]);       % no. of nodes 
        n_dim = uint32([]);         % no. of spatial dimensions
        n_dof = uint32([]);         % no. of global degrees of freedom
        n_dof_node = uint32(1);     % no. of dofs per node (scalar = 1)  
        nodes = [];                 % global coordinates (no. nodes, no. dim)
        node_plot_order;            % node plotting order (only needed if nodes are out of order)
        opt = ...                   % Options structure
          struct('space', 'scalar');
    end
    
    % Public properties (read only; except FEmesh and Element)
    properties (GetAccess = public, SetAccess = {?mFEM.FEmesh, ?mFEM.Element})
        on_boundary;                % flag if element is on a boundary
        boundary_id = uint32([]);   % list of all boundary ids for element
        subdomain = uint32([]);     % list of all subdomain ids for element
        side;                       % side info, see constructor        
        local_n_dim;                % local dimensions
        direct = false;             % a flag for using direct assembly, see Truss
    end
    
    % Protected properties
     properties (Hidden = true, Access = {?mFEM.FEmesh, ?mFEM.Element}) 
         neighbors;       % storage of nieghbor elements (see FEmesh.find_neighbors)
         global_dof = []; % Global dof for nodes of element
                          % (these are not the true dofs (except in scalar 
                          % space) as such the user should always access the
                          % dofs for an element with the get_dof() function.) 
     end    
   
    % Public Methods
    % (These methods are accessible by the user to create the element and
    % access the shape functions and other necessary parameters)
    methods (Access = public)
        function obj = Element(id, nodes, varargin)
            %ELEMENT Class constructor.
            %
            % This is an abstract class, it must be inherited by a subclass
            % to operate, see Line2.m for example. The following syntax and
            % descriptions apply to all subclasses unless noted otherwise
            % in the documentation for the specific element.
            %
            % Syntax
            %   Element(id, nodes)
            %   Element(id, nodes, 'PropertyName', PropertyValue, ...)
            %
            % Description
            %   Element(id, nodes) creates an element given, where id is a
            %   unique identification number for this element and nodes is 
            %   a matrix of node coordinates (global) that should be 
            %   arranged as column matrix (no. nodes x no. dims).
            %
            %   Element(id, nodes, 'PropertyName', PropertyValue, ...) 
            %   allows the user to customize the behavior of the element, 
            %   the available properties are listed below.
            %
            % Element Property Descriptions
            %   space
            %       {'scalar'} | 'vector'  | integer
            %       Allows the type of FEM space to be set: scalar sets the 
            %       number of dofs per node to 1, vector  sets it to the 
            %       no. of space dimension, and  specifing a number sets it
            %       to that value.

            % Insert required values into object properties
            obj.id = id;
            obj.nodes = nodes;
            [obj.n_nodes, obj.n_dim] = size(nodes);
            obj.local_n_dim = obj.n_dim;

            % Change dofs per node
            if nargin == 4 && strcmpi(varargin{2},'vector');
                obj.n_dof_node = obj.n_dim;  
                obj.opt.space = varargin{2};
            elseif nargin == 4 && isnumeric(varargin{2});
            	obj.n_dof_node = varargin{2};
                obj.opt.space = varargin{2};
            end   
            
            % Determine the total number of global dofs
            obj.n_dof = obj.n_nodes * obj.n_dof_node;
            
            % Intialize the side data structure
            obj.side = struct('on_boundary', true, ...
                'boundary_id', cell(obj.n_sides,1),...
                'neighbor',[], 'neighbor_side', []);
            
            % Initialize neighbor array
            obj.neighbors = feval([class(obj),'.empty']);
        end
        
        function N = shape(obj, varargin)
            %SHAPE Returns the shape functions
            %
            % Syntax
            %   shape(xi)
            %   shape(xi,eta)
            %   shape(xi,eta,zeta)
            %   shape(...,'PropertyName',PropertyValue)
            %
            % Description
            %   shape(...) returns the element shape functions evaluated at
            %   the locations specified in the inputs, the number of which
            %   varies with the number of space dimensions.
            %
            %   shape(...,'PropertyName',PropertyValue) allows user to
            %   override the vectorized output using the scalar flag, this
            %   is used by GET_POSITION
            %
            % SHAPE Property Descriptions
            %   scalar
            %       true | {false}
            %       If this is set to true the vectorized output is
            %       ignored. E.g.,
            %           N = shape(0,'-scalar');

            % Parse options (do not use gather_user_options for speed)
            options.scalar = false;
            if strcmpi(varargin{end},'-scalar');
                varargin = varargin(1:end-1);
                options.scalar = true;
            elseif nargin > 2 && ischar(varargin{end-1});
                options.scalar = varargin{end};
                varargin = varargin(1:end-2);             
            end                
            
            % Scalar field basis functions
            N = obj.basis(varargin{:});

            % Non-scalar fields
            if ~options.scalar && (obj.n_dof_node > 1 && strcmpi(obj.opt.space,'vector'));
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
            %SHAPE_DERIV Returns shape function derivatives in global x,y system
            %
            % Syntax
            %   shape_deriv(xi)
            %   shape_deriv(xi,eta)
            %   shape_deriv(xi,eta,zeta)
            %
            % Description
            %   shape_deriv(...) returns the element shape function 
            %   derivatives evaluated at the locations specified in the 
            %   inputs, the number of which varies with the number of space 
            %   dimensions.

            % Scalar field basis functin derivatives
            B = obj.grad_basis(varargin{:});
                        
            % Non-scalar fields
            if obj.n_dof_node > 1 && strcmpi(obj.opt.space,'vector');
                b = B;                      % Re-assign scalar basis
                r = obj.n_dof_node;         % no. of rows
                c = r*size(b,2);            % no. of cols
                B = zeros(r+1,c);           % size the vector basis

                % Loop through the rows and assign scalar basis
                for i = 1:r;
                    B(i,i:r:c)  = b(i,:);
                    B(r+1, i:r:c) = b((r+1)-i,:);
                end
           
%             % Truss space    
%             elseif strcmpi(obj.opt.space, 'truss');
%                 b = B;
%                 B = zeros(1,obj.n_nodes*obj.n_dim);
%                 B(1:obj.n_dim:end) = b; 
            end
        end
        
        function size(obj)
            %SIZE Returns the length, area, or volume of the element
            %
            % Syntax
            %   L = size()
            %
            % Description 
            %   L = size() returns the size of the element, the meaning of
            %   the size depends on the element, 1D, 2D, or 3D, which
            %   should return the length, area, or volume, respectively.
            %
            % This function must be defined within the subclass, e.g., see
            % Line2. If it is not defined you will recieve an error if an
            % attempt use is made. This was not made an abstract method,
            % because in many cases this property is not needed. However,
            % the System class has a variable reserved L for it, so it
            % should be available for all Elements, even if it is an error.
            %
            % See Also
            %   Line2
            
            % The default behavior is an error message
            error('Element:size','The size method is not defined for the %s element type',class(obj));
        end
        
        function L = hmax(obj)
            %HMAX Return maximum length between ANY two points
            %
            % Syntax
            %   L = max();
            %
            % Description
            %   L = max() returns the maximum distance between any two
            %   points on the element.
            L = max(obj.distance());
        end
        
        function L = hmin(obj)
            %HMIN Return inimum length between ANY two points
            %
            % Syntax
            %   L = min();
            %
            % Description
            %   L = min() returns the minimum distance between any two
            %   points on the element.
            L = min(obj.distance());
        end
            
        function J = detJ(obj, varargin)
            %DETJ Returns the determinate of the jacobian matrix
            %
            % Syntax
            %   detJ(xi)
            %   detJ(xi,eta)
            %   detJ(xi,eta,zeta)
            %
            % Description
            %   detJ(...) returns the determinante of the jacobian 
            %   evaluated at the locations specified in the inputs, the 
            %   number of which varies with the number of space dimensions.
            J = det(obj.jacobian(varargin{:}));
        end 
        
        function varargout = get_position(obj, varargin)
            %GET_POSITION Returns the real coordinates given xi, eta, ...
            %
            % Syntax
            %   x = get_position(xi)
            %   [x,y] = get_position(xi,eta)
            %   [x,y,z] = get_position(xi,eta,zeta)
            %   [...] = get_position(..., 'PropertyName', PropertyValue)
            %   xyz = get_position(...)
            %
            % Description
            %   [...] = get_position(...) returns the position in global
            %   system given a value(s) for the local position, the number 
            %   of outputs varies according the number of spatial 
            %   dimensions.
            %
            %   [...] = get_position(..., 'PropertyName', PropertyValue)
            %   allow the limiting of which shape functions are used for
            %   the mapping, see EXAMPLE14b for an example.
            %
            %   xyz = get_position(...) same as above but it returns the
            %   positions as a single array.
            %
            % GET_POSITION Property Descriptions
            %   index
            %       logical array | numeric array
            %       An array that limits what shape functions are used for
            %       the mapping, this is useful elements that have mutiple
            %       dofs per node. The value of index should be defined
            %       such that the following is valid
            %           N = obj.shape(xi,...)
            %           N = N(opt.index)
            %       
            %       Using this property overrides the default behavior,
            %       which is to use opt.index = 1:n_dof_node:end when a the
            %       no. of nodes is different from no. of shape functions,
            %       as is the case for the Beam element. For the Beam
            %       element there are 2 dofs per node and the 1 and 3 value
            %       for the shape functions may be used for maping the
            %       displacement, this behavior is the default. If you have
            %       an element that behaves otherwise, use the index
            %       property.
            
            if strcmpi(obj.opt.space,'vector');
                warning('Element:get_position','This is not tested for vector space and likely does not work correctly');
            end
            
            % Set the default options
            options.index = [];
            
            % Parse inputs
            for i = 1:length(varargin);
                TF = ischar(varargin{i});
                if TF;
                    options = gather_user_options(options,varargin{i:end});
                    varargin = varargin(1:i-1);
                    break;  
                end  
            end

            % The number of spatial dimensions available
            n = obj.n_dim;

            % Loop through the dimensions and return the desired position
            xyz = zeros(1,n);
            for i = 1:n
                
               % Evaluate shape functions
               N = obj.shape(varargin{:},'-scalar');

               % The index property was set, so limit the shape functions
               if ~isempty(options.index);
                   N = N(options.index);
                   
               % The no. of shape functions is not equal to no. of nodes,
               % so assume that the 1:n_dof_node:end shape functions are used to
               % map the position, this is the case for Beam elements.
               elseif length(N) > length(varargin);
                   N  = N(1:obj.n_dof_node:end);
               end

               % Map the position
               xyz(i) = N*obj.nodes(:,i);
            end
   
            % Reduce to single array if only a single output is given
            if nargout > 1;
                varargout = num2cell(xyz,1);
            else
                varargout{1} = xyz;
            end
        end
        
        function n = get_normal(obj, varargin)
            %GET_NORMAL Returns the normal vector at a location
            %
            % Syntax
            %   n = get_normal(xi)
            %   n = get_normal(xi,eta)
            %
            % Description
            %   n = get_normal(...) returns the normal vector for a given
            %   local location. The number of inputs depends on the number 
            %   of space dimensions of the element. 
            %  
            % This function is still under development and is not fully
            % tested, expect changes.
            
            % The number of spatial dimensions available
            n = obj.n_dim; 
            
            % 2D: Side is defined by a line
            if n == 2;
                % Compute the tangent at the point
                n = obj.local_grad_basis(varargin{:}) * obj.nodes;
                
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
            %BUILD_SIDE Build an element for the side
            %
            % Syntax
            %   side = build_side(id)
            % 
            % Description
            %   side = build_side(id) creates an element for the specified
            %   side, the type of element is specified in the side_type
            %   property.
            
            if obj.n_dim == 3;
                warning('Element:build_side','Feature not tested in 3D');
            end
            
            % Extract the nodes for the side
            dof = obj.side_dof(id,:);
            node = obj.nodes(dof,:);

            % Create the side element
            side = feval(['mFEM.elements.',obj.side_type], NaN, node, ...
                'Space', obj.opt.space);
            
            % Set the global dofs
            side.global_dof = obj.global_dof(dof);
        end
               
        function dof = get_dof(obj, varargin)
            %GET_DOF The global degrees of freedom, account for type of space
            %
            % Syntax
            %   dof = get_dof()
            %   dof = get_dof('PropertyName', PropertyValue,...)
            % 
            % Description           
            %   dof = get_dof() returns all of the global dofs for element
            %
            %   dof = get_dof('PropertyName', PropertyValue) returns the
            %   dofs for the element subject to the restrictions specified
            %   by the properties, see the descriptions below for details.           
            %
            % GET_DOF Property Descriptions
            %   Side
            %       integer
            %       Indicates that only the degrees of freedom for the
            %       specific side should be returned.
            %
            %   Local
            %       true | {false}
            %       Toggles the type of degrees of freedom to return, which
            %       is generally import for sides. For example, the
            %       following returns the local degrees of freedom for side
            %       number 1 of an element.
            %           dof = get_dof('Side',1,'local',true) or
            %           dof = get_dof('Side',1,'-local')
            
            % Set default options and gather the options
            options.side = [];
            options.local = false;
            options = gather_user_options(options, varargin{:});
            
            % Extract ALL of the degrees of freedom
            if isempty(options.side)
                if options.local
                    dof = 1:obj.n_dof;
                else
                    dof = obj.global_dof;
                end
            
            % Extract dofs for the specified side    
            elseif isnumeric(options.side) && options.side <= obj.n_sides;
                s = options.side;
                if options.local;
                    dof = obj.side_dof(s,:);
                else
                    idx = obj.side_dof(s,:);
                    dof = obj.global_dof(idx);
                end
            else
                error('Element:get_dof','Input error.');
            end

            % Non-scalar FE space
            if obj.n_dof_node > 1;
                dof = obj.transform_dof(dof);
            end
        end  
    end
    
    methods (Hidden = true, Access = private)        
        function D = transform_dof(obj, d)
            %TRANSFROM_DOF Converts the dofs for vector element space
            %
            % Syntax
            %   D = transform_dof(d)
            %   
            % Description
            %   D = transform_dof(d) converts the scalar degrees of freedom
            %   for to vector based degrees of freedom. For example,
            %   inputing d = [1,3] returns D = [1,2,5,6].
            
            D = transform_dof(d,obj.n_dof_node); % call func in bin
        end
        
        function d = distance(obj)
            %DISTANCE Compute distances between all points
            %
            % Syntax
            %   d = distance()
            %
            % Description
            %   d = distance() computes the distance between all nodes, so
            %   is returns N! distances, where N is the number of nodes
            %   minus one.
           
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
    