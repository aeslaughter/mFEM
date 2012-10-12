classdef Linear2 < mFEM.Element
    %Linear 2-node 1D linear element
    %
    %      (-1)   (1)   (1)
    %         1---------2
    %

    % Define the inherited abstract properties
    properties (SetAccess = protected, GetAccess = public)
        n_shape = 2;                % no. of shape functions
        n_sides = 2;                % no. "sides" (nodes are sides in 1D)
        side_dof = [1; 2];          % local dofs of the "sides"
        side_defn = [1,-1; 1,1];    % xi definitions for "sides"
    end
    
    % Define the Linear2 constructor
    methods 
        function obj = Linear2(id, nodes, varargin)
           % Class constructor; calls base class constructor
           
           % Test that nodes is sized correctly
           if ~all(size(nodes) == [2,1]);
                error('Linear2:Linear2','Nodes not specified correctly; expected a [2 x 1] array, but recieved a [%d x %d] array.', size(nodes,1), size(nodes,2));
           end
           
           % Call the base class constructor
           obj = obj@mFEM.Element(id, nodes, varargin{:}); 
        end
    end
    
    % Define the inherited abstract methods (protected)
    methods (Access = protected)      
        function N = basis(~, xi)
            % Returns a row vector of local shape functions
            N(1) = 1/2*(1-xi);
            N(2) = 1/2*(1+xi);
        end

        function B = grad_basis(obj, varargin) 
            % Gradient of shape functions
            
            % grad(B) is constant, produce a warning for agruements
            if nargin > 1;
                warning('Tri3:grad_basis', 'The shape function deriviatives are constant for the Linear2 element, thus no spatial coordinates are needed.');
            end
            
            % Proper gradient
            B = inv(obj.jacobian()) * [-1/2, 1/2];
        end
             
        function J = jacobian(obj, varargin)
            % Returns the jacobian matrix  
            
            % The Jacobian is constant, produce a warning for agruements
            if nargin > 1;
                warning('Tri3:jacobian', 'The Jacobian for the Linear2 element is constant, thus no spatial coordinates are needed.');
            end

            % Return the Jacobian matrix
            J = [-1/2, 1/2]*obj.nodes;                 
        end
    end
end