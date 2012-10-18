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
      n_sides;      % no. of sides
      lims;         % limits of local coordinate; assumesall dim. vary the same
      side_dof;     % array of local node ids for each side
      side_type;    % defines the type of element that defines the sides.
    end
    
    % Abstract Methods (protected)
    % (the user must redfine these in subclasse, e.g. Quad4)
    methods (Abstract, Access = protected)
        N = basis(obj, varargin)          % basis functions
        B = grad_basis(obj, varargin)     % basis function derivatives (dN/dx, ...)
        G = local_grad_basis(obj, varargin) % basis function derivatives (dN/dxi, ...)
        J = jacobian(obj, varargin)       % the Jacobian matrix for the element
    end
    
    % Public properties (read only)
    properties (SetAccess = protected, GetAccess = public)
        id = [];          % element id [double]
        nodes = [];       % global coordinates (no. nodes, no. dim) [double]
        n_nodes = [];     % no. of nodes [double]
        n_dim = []        % no. of spatial dimensions [double]
        n_dof = [];       % no. of global degrees of freedom
        space = 'scalar'; % scalar or vector space [string]
        is_side = false;  % id's the element as a side or not
    end
    
    % Public properties (read only; except FEmesh)
    properties (SetAccess = {?mFEM.FEmesh, ?mFEM.Element}, SetAccess = protected, GetAccess = public)
        on_boundary;                % flag if element is on a boundary
        boundary_id = uint32([]);   % list of all boundary ids for element
        neighbors;                  % list of all neighbor elements (will replace side.neighbor)
        n_neighbors = [];           % no. of neighbor elements

        % structure containing side info
        side = struct('on_boundary', [], 'boundary_id', uint32([]),...
            'dof', uint32([]), 'global_dof', uint32([]), ...
            'neighbor', uint32([]));          
    end
    
    % Protected properties
    properties (Access = {?mFEM.FEmesh, ?mFEM.Element}, Access = protected)
       global_dof = []; % global dof for nodes of element    
       side_nodes = []; % nodal coord. for side elements, see get_normal
%       n_dof_node = 1;  % no. of dofs per node 
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
            %   Element(id, nodes, space)
            %
            % Creates an element given:
            %   id: unique identification number for this element
            %   nodes: matrix of node coordinates (global), it should be 
            %          arranged as a column matrix [no. nodes x no. dims]
            %   space (optional): string that is 'scalar' (default) or
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
                obj.n_dof_node = obj.n_dim;
            end
            
%             % Initialize side data structure
%            	n = obj.n_nodes;
%             obj.side(n).on_boundary = [];
%             obj.side(n).boundary_id = uint32([]);
%          %   obj.side(n).has_neighbor = [];
%             obj.side(n).dof = uint32([]);
%             obj.side(n).global_dof = uint32([]);
%             obj.side(n).neighbor = uint32([]); %struct('element', uint32([]), 'side', uint32([]));;
        end
        
        function N = shape(obj, varargin)
            % Returns the shape functions

            % Scalar field basis functions
            N = obj.basis(varargin{:});

            % Vector
            if strcmpi(obj.space, 'vector');
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
           
            % Initialize the output
            varargout = cell(obj.n_dim);
            
            % Determine the nodes based on if the element is a side or not
            if obj.is_side;
                n = obj.n_dim + 1;
                node = obj.side_nodes;
            else
                n = obj.n_dim;
                node = obj.nodes;
            end

            % Loop through the dimensions and return the desired position
            for i = 1:n;
               varargout{i} = obj.shape(varargin{:})*node(:,i); 
            end

        end
        
        function n = get_normal(obj, varargin)
            % Returns the normal vector at a given xi, eta, ...
            %
            % The element must be a side element created with build_side
            % of a parent element for this function to be used.
           
            % Throw an error if the element is not a side
            if ~obj.is_side;
                error('Element:get_normal', 'Function only available for side elements');
            end
            
            % 1D: Side is defined by a line
            if obj.n_dim == 1;
                % Compute the tangent at the point
                n = obj.local_grad_basis(varargin{:}) * obj.side_nodes;
                
                % Re-arrange tangent to give the normal (outward from
                % element face is positive)
                n = [n(2), -n(1)]/norm(n);

            % 2D: Side is defined by a plane    
            elseif obj.n_dim == 2;
                error('Element:get_normal', 'Not yet supported');
                
            % Only defined for 1D and 2D sides
            else
                error('Element:get_normal', 'Not defined for %d-D side element', obj.n_dim);
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
            
            % Create a map to the line/plane
            mapped = zeros(size(node,2),1);
            for i = 2:size(node,1);
                mapped(i,:) = norm(node(i,:) - node(1,:));
            end
            
            % Create the side element and flag it as a side
            side = feval(['mFEM.',obj.side_type], NaN, mapped, obj.space);
            
            % Indicate that the element created was a side
            side.is_side = true;
            
            % Store the actual coordinates in the side element
            side.side_nodes = node;
            
            % Override the no. of dofs per node to match the parent element
 %           side.n_dof_node = obj.n_dof_node;
        end
               
        function dof = get_dof(obj)
            % The global degrees of freedom, account for type of space
            
            % Scalar FE space
            if strcmpi(obj.space,'scalar');
                dof = obj.global_dof;
                
            % Vector FE space
            elseif strcmpi(obj.space,'vector');
                dof = obj.transform_dof(obj.global_dof);
            else
                error('Element:get_dof', 'Unknown finite elment space, %s, for element %d\n', obj.space, obj.id); 
            end
        end  
        
        function dof = get_side_dof(obj, s)
            % Extract the local dofs for the specified side
           
            % Scalar FE space
            if strcmpi(obj.space,'scalar');
                dof = obj.side(s).dof;

            % Vector FE space
            elseif strcmpi(obj.space,'vector');
                dof = obj.transform_dof(obj.side(s).dof);
            else
                error('Element:get_dof', 'Unknown finite elment space, %s, for element %d\n', obj.space, obj.id); 
            end
        end
        
        function [s_e, s_n] = match_sides(elem, neighbor)
            % Returns the sides that match an element and its neighbor
            
            % Check that neighbor is a neighbor 
            if ~any(neighbor == elem.neighbors);
               error('Element:match_sides', 'The neighbor element (%d) supplied is not a neighbor of the current element (%d)', neighbor.id, elem.id); 
            end
            
            % Extract the global dofs from the element
            a = elem.nodes;
            b = neighbor.nodes;
           
            % Locate where a & b match
            [~,ia,ib] = intersect(a,b,'rows','R2012a');
  
            % Locate the matching side on the current element
            C = a(sort(ia),:);
            for s = 1:elem.n_sides;
                if isequal(C, elem.nodes(sort(elem.side_dof(s,:)),:));
                   s_e = s;
                   break;
                end
            end
            
            % Locate the matching side on the neighbor element
            C = b(sort(ib),:);             
            for s = 1:neighbor.n_sides;
                if isequal(C, neighbor.nodes(sort(neighbor.side_dof(s,:)),:));
                   s_n = s;
                   break;
                end
            end  
        end
        
        function D = transform_dof(obj, d)
            % Converts the dofs for vector element space
            
            n = obj.n_dim;              % no. of dimensions
            D = zeros(n*length(d),1);   % size of vector space dofs
            
            % Loop through dimensions and build vector
            for i = 1:n;
                D(i:n:end) = d*n - (n-i);
            end 
        end
    end
end
    