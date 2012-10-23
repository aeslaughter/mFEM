classdef Truss2 < mFEM.Element
    % A 2-node 1D linear element, but located in 2D space
    %
    %      (-1)   (1)   (1)
    %         1---------2
    %

    % Define the inherited abstract properties
    properties (SetAccess = protected, GetAccess = public)
        n_sides = 2;                % no. "sides" (nodes are sides in 1D)
        lims = [-1,1];              % limits xi
        side_dof = [1; 2];          % local dofs of the "sides"
        side_type = '';             % 1D elements do not have side elements
    end
    
    % Define the Linear2 constructor
    methods 
        function obj = Truss2(id, nodes, varargin)
           % Class constructor; calls base class constructor
           
           % Test that nodes is sized correctly
           if ~all(size(nodes) == [2,2]);
                error('Truss2:Truss2','Nodes not specified correctly; expected a [2x2] array, but recieved a [%dx%d] array.', size(nodes,1), size(nodes,2));
           end
           
           % Call the base class constructor
           obj = obj@mFEM.Element(id, nodes, varargin{:}); 
           
           % Indicate that the local dimension is different from spatial
           obj.local_dim = 1;
        end
        
        % Define the size function
        function L = size(obj)
        	L = norm(obj.nodes);   
        end
    end
    
    % Define the inherited abstract methods (protected)
    methods (Access = protected)      
        function N = basis(~, xi)
            % Returns a row vector of local shape functions
            N(1) = 1/2*(1 - xi);
            N(2) = 1/2*(1 + xi);
        end

        function B = grad_basis(obj, varargin) 
            % Gradient of shape functions
            
            % grad(B) is constant, produce a warning for agruements
            if nargin > 1;
                warning('Truss2:grad_basis', 'The shape function deriviatives are constant for the Linear2 element, thus no spatial coordinates are needed.');
            end
            
            % Proper gradient
            B = inv(obj.jacobian()) * obj.local_grad_basis;
        end
             
        function J = jacobian(obj, varargin)
            % Returns the jacobian matrix  
            
            % The Jacobian is constant, produce a warning for agruements
            if nargin > 1;
                warning('Truss2:jacobian', 'The Jacobian for the Truss2 element is constant, thus no spatial coordinates are needed.');
            end

            % Return the Jacobian matrix (1/2 the length)
            J = 1/2 * obj.size();               
        end
        
        function G = local_grad_basis(obj, varargin)
            G = [-1/2, 1/2];
        end
    end
end