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
    
    % Abstract Properties
    properties (Abstract = true, SetAccess = protected, GetAccess = public)  
      n_shape;     % no. of shape functions
      n_sides;     % no. of sides
      sides;       % array of xi, eta values for sides
      side_nodes;  % array of local node ids for each side
    end
    
    % Abstract Methods (protected)
    % (the user must redfine this in subclasses, e.g. Quad4)
    methods (Abstract, Access = protected)
         N = basis(~, varargin)          % local basis functions
        GN = grad_basis(~, varargin)     % local basis function derivatives
    end
    
    % Public properties (read only)
    properties (SetAccess = protected, GetAccess = public)
        id = [];          % element id   
        nodes = [];       % global coordinates [no. nodes, no. dim]
        n_nodes = [];     % no. of nodes
        n_dim = []        % no. of spatial dimensions
        space = 'scalar'; % scalar or vector space
    end
    
    % Public properties (read only; except FEmesh)
    properties (SetAccess = {?FEmesh}, SetAccess = protected, GetAccess = public)
        % structure containing neighbor information
        neighbor = struct('element',[],'side_index',[]);                
    end
   
     % Private properties (except FEmesh)
     properties (SetAccess = {?FEmesh}, SetAccess = protected, GetAccess = public)
     	global_dof = []; % vector of global dof for the nodes of this element
     end
    
    % Public Methods
    % These methods are accessible by the user to create the element and
    % access the shape functions and other necessary parameters
    methods
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
        end
        
        function N = shape(obj, varargin)
            % Returns the shape functions
            
            % Scalar field basis functions
            N = obj.basis(varargin{:});
            
            % Vector
            if strcmpi(obj.type, 'vector');
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
            B = inv(obj.jacobian(xi,eta)) * obj.grad_basis(varargin{:});
            
            % Vector
            if strcmpi(obj.type, 'vector');
                b = B;                      % Re-assign scalar basis
                r = obj.n_dim;              % no. of rows
                c = obj.n_dim*length(B);    % no. of cols
                B = zeros(r+1,c);           % size the vector basis
                
                % Loop through the rows and assign scalar basis
                for i = 1:r;
                    B(i,i:r:c) = b(i,:);
                    B(r+1, i:r:c) = b(i,:);
                end
            end
        end
        
        function N = side_shape(obj, id, varargin)
            % Returns the shape functions the side identified by id
            side = obj.sides(id, :);
            side(isnan(side)) = varargin{:};
            side = num2cell(side);
            N = obj.shape(side{:});
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
                error('ERROR: Unknown finite elment space (%s) for element %d\n', obj.space, obj.id);
                
            end
        end
        
        function J = detJ(obj, varargin)
            % Returns the determinate of the jacobian matrix
            J = det(obj.jacobian(varargin{:}));
        end 
      
        function J = jacobian(obj, varargin)
            % Returns the jacobian matrix
            J = obj.grad_basis(varargin{:})*obj.nodes;                    
        end       
    end
end
    