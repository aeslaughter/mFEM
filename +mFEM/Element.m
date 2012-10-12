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
        N = basis(obj, varargin)          % local basis functions
        B = grad_basis(obj, varargin)     % local basis function derivatives
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
    end
    
    % Public properties (read only; except FEmesh)
    properties (SetAccess = {?mFEM.FEmesh, ?mFEM.Element}, SetAccess = protected, GetAccess = public)
        side = struct([]);          % structure containing side info
        on_boundary;                % flag if element is on a boundary
        boundary_id = uint32([]);   % list of all boundary ids for element
      	global_dof = [];            % global dof for nodes of element       
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
            end
            
            % Initialize side data structure
           	n = obj.n_nodes;
            obj.side(n).boundary_id = uint32([]);
            obj.side(n).has_neighbor = [];
            obj.side(n).dof = uint32([]);
            obj.side(n).global_dof = uint32([]);
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
        
        function side = build_side(obj, id)
            % Build an element for the side
            
            if obj.n_dim == 3;
                error('Element:build_side','Feature not yet supported');
            end
            
            % Extract the nodes for the side
            dof = obj.side_dof(id,:);
            node = obj.nodes(dof,:);
            
            % Create a map to the line/plane
            mapped = zeros(size(node),1);
            for i = 2:size(node,1);
                mapped(i,:) = norm(node(i,:) - node(1,:));
            end
            
            % Create the side element and flag it as a side
            side = feval(['mFEM.',obj.side_type], NaN, mapped, obj.space);
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
end
    