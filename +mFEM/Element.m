classdef Element < handle
    %Element Base class for defining elements
    % Inludes the general behavior of an element, including the node locations,
    % id, shape functions, etc...
    %
    % This is an abstract class, as such it must be inherited to function.
    % The abstract properties and methods must be redifined in the
    % subclass, see Quad4.m for an example.
    %
    % see Quad4
    
    % Abstract Properties (must be redefined in subclass)
    properties (Abstract = true, SetAccess = protected, GetAccess = public)  
      n_shape;     % no. of shape functions
      n_sides;     % no. of sides
      
      % Array of local node ids for each side, should be an m x n array,
      % where m = number of sides and n = number of nodes on the side. See
      % Quad4.m for an example
      side_dof;    % array of local node ids for each side
      
      % An array defining the values of the local coordinate system that
      % are fixed for the sides. It should be an m x 2 array, where m is
      % the number of sides. The first column gives the index of the local
      % coordinate system (e.g., 1 = xi and 2 = eta; see Quad4.m). The
      % second column defines the value that the coordinate is fixed. For
      % example, a row of [1,-1] indiates that the side is defined when xi
      % = -1 (see Quad4.m)
      side_defn;   % array defining the values of xi or eta that are fixed
      
      % It is possible to create an element with nodes that are the same,
      % as in the case when triangle element is generated from a
      % quadraleral. The collapsed_nodes array indicates which nodes are
      % the same. See Tri3 for an example. If all nodes are unqiue make
      % this an empty array.
      collapsed_nodes;   % lists the nodes that are collapsed (see Tri3)
    end
    
    % Abstract Methods (protected)
    % (the user must redfine this in subclasses, e.g. Quad4)
    methods (Abstract, Access = protected)
         N = basis(~, varargin)          % local basis functions
        GN = grad_basis(~, varargin)     % local basis function derivatives
    end
    
    % Public properties (read only)
    properties (SetAccess = protected, GetAccess = public)
        id = [];          % element id [double]
        nodes = [];       % global coordinates (no. nodes, no. dim) [double]
        n_nodes = [];     % no. of nodes [double]
        n_dim = []        % no. of spatial dimensions [double]
        n_dof = [];       % no. of global degrees of freedom
        space = 'scalar'; % scalar or vector space [string]
    end
    
    % Public properties (read only; except FEmesh)
    properties (SetAccess = {?mFEM.FEmesh, ?Element}, SetAccess = protected, GetAccess = public)
        side = struct([]); % structure containing side information (computed by FEmesh)
        on_boundary; % flag if element is on a boundary (has a side w/o neighbor) [bool]
        boundary_id = uint32([]); % list of all boundary ids for the element
    end
   
    % Private properties (except FEmesh)
    properties (Access = {?mFEM.FEmesh, ?Element}, Access = protected)
     	global_dof = [];        % vector of global dof for the nodes of this element
     end
    
    % Public Methods
    % These methods are accessible by the user to create the element and
    % access the shape functions and other necessary parameters
    methods (Access = public)
        function obj = Element(id, nodes, varargin)
            % Class constructor.
            %
            % This is an abstract class, it must be inherited by a subclass
            % to operate, see Quad4.m for example.
            %
            % Syntax:
            %   Element(id, nodes)
            %   Element(id, nodes, type)
            %
            % Creates an element given:
            %   id: unique identification number for this element
            %   nodes: matrix of node coordinates (global), it should be 
            %          arranged as a column matrix [no. nodes x no. dims]
            %   type (optional): string that is 'scalar' (default) or
            %       'vector', which indicates the type of solution
            %
            
            % Insert required values into object properties
            obj.id = id;
            obj.nodes = nodes;
            [obj.n_nodes, obj.n_dim] = size(nodes);
    
            % User supplied the space
            if nargin == 3;
                obj.space = varargin{1};
            end         
            
            % Determine the number of global dofs
            obj.n_dof = obj.n_nodes;
            if strcmpi(obj.space,'vector');
                obj.n_dof = obj.n_dof * obj.n_dim;
            end
            
            % Initialize side data structure
           	n = obj.n_nodes;
            obj.side(n).boundary_id = uint32([]);
            obj.side(n).global_dof = uint32([]);
            obj.side(n).has_neighbor = [];
            obj.side(n).neighbor = struct('element', uint32([]), 'side', uint32([]));;
        end
        
        function N = shape(obj, varargin)
            % Returns the shape functions
            
            % Scalar field basis functions
            N = obj.basis(varargin{:});
            
            % Vector
            if strcmpi(obj.space, 'vector');
                n = N;                      % Re-assign scalar basis
                r = obj.n_dim;              % no. of rows
                c = obj.n_dim*length(N);    % no. of cols
                N = zeros(r,c);             % size the vector basis
                
                % Loop through the rows and assign scalar basis
                for i = 1:r;
                    N(i,i:r:c) = n;
                end
            end            
        end
        
        function B = shape_deriv(obj, varargin)
            % Returns the shape function derivatives in x,y system

            % Scalar field basis functin derivatives
            B = inv(obj.jacobian(varargin{:})) * obj.grad_basis(varargin{:});
            
            % Adjust for collapsed nodes
            if ~isempty(obj.collapsed_nodes);
                idx = true(size(B,2),1);
                c = obj.collapsed_nodes(:,1);
                idx(c) = false;
                B = B(:,idx);
            end
             
            % Vector
            if strcmpi(obj.space, 'vector');
                b = B;                      % Re-assign scalar basis
                r = obj.n_dim;              % no. of rows
                c = r*size(b,2);            % no. of cols
                B = zeros(r+1,c);           % size the vector basis

                % Loop through the rows and assign scalar basis
                for i = 1:r;
                    B(i,i:r:c)  = b(i,:);
                    B(r+1, i:r:c) = b((r+1)-i,:);
                end
            end
        end
        
        function J = detJ(obj, varargin)
            % Returns the determinate of the jacobian matrix
            J = det(obj.jacobian(varargin{:}));
        end 
        
        function N = side_shape(obj, id, varargin)
            % Returns the shape functions the side identified by id
            in = obj.side_input(id, varargin{:});
            N = obj.shape(in{:});
        end
        
        function B = side_shape_deriv(obj, id, varargin)
           % Returns the shape function derivatives for the specified side
           in = obj.side_input(id, varargin{:});
           B = obj.shape_deriv(in{:});
        end
        
        function J = side_detJ(obj, id, varargin)
            % Return derterminant of Jacobian for the specified side
            
            % Collect input in proper form
            in = obj.side_input(id, varargin{:});
            
            % Get the nodes for the system
            nodes = obj.get_nodes();
            
            % Local dofs for the current side
            idx = obj.side_dof(id,:);
          
            % Perform calculation based on dimension of system
            if obj.n_dim == 1;
                error('Element:side_detJ', 'Not defined for 1D elements.');

            % 2D: The side is a line, so detJ = length/2
            elseif obj.n_dim == 2;
                % Use max distance between any of the nodes
                J = max(pdist(nodes(idx,:)))/2;
                
                if ~strcmpi('mFEM.Quad4', class(obj));
                    warning('Element:side_detJ', 'The method for computing the length/2 (Jacobian of line) was only tested with a Quad4 element.');
                end
                
            % 3D: The side is a face    
            else
                % Compute Jacobain for the side
                GN = obj.grad_basis(in{:});
                J = det(GN(:,idx)*nodes(idx,:));
                
                warning('Element:side_detJ', 'The method for computing the Jacobian on sides of 3D elements was not tested.');
            end
        end
        
        function dof = get_dof(obj)
            % The global degrees of freedom, account for type of space
            
            % Scalar FE space
            if strcmpi(obj.space,'scalar');
                dof = obj.global_dof;
                
            % Vector FE space
            elseif strcmpi(obj.space,'vector');
                D = obj.global_dof;
                dof = [];
                for i = 1:length(D);
                    d0 = 2*(D(i) - 1) + 1;
                    dn = d0 + obj.n_dim-1;
                    dof = [dof; (d0:dn)'];
                end
                
            else
                error('Element:get_dof', 'Unknown finite elment space (%s) for element %d\n', obj.space, obj.id); 
            end
        end  
    end
    
    % Private methods
    methods (Access = private)
        
        function nodes = get_nodes(obj)
            % Retrn a matrix of the nodes, accounting for collapsed nodes
            
            % The unique nodes for the element
            nodes = obj.nodes;
            
            
            % Account for  collapsed nodes
            if ~isempty(obj.collapsed_nodes);
                ix1 = obj.collapsed_nodes(:,1);
                ix2 = obj.collapsed_nodes(:,2);
                nodes(ix1,:) = obj.nodes(ix2,:);
            end  
        end
        
        function J = jacobian(obj, varargin)
            % Returns the jacobian matrix  
            
            J = obj.grad_basis(varargin{:})*obj.get_nodes();                    
        end
        
        function in = side_input(obj, id, varargin)
            % Parses side input for use in shape and shape_deriv
            
            % Extract index and value for current side
            s = obj.side_defn(id, :);
            
            % 1D case, there is no input
            if obj.n_dim == 1;
               in = num2cell(s(2)); 
               return;
            end
            
            % Create correctly sized index vector and input vector 
            idx = 1:obj.n_dim;
            in(idx) = 0;
            
            % Insert correct values 
            in(idx ~= s(1)) = varargin{:};
            in(s(1)) = s(2);
 
            % Convert numeric array to cell array
            in = num2cell(in);
        end
    end
end
    