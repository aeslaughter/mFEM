classdef Line3 < mFEM.Element
    %LINE3 A 3-node, 1D quadradtic element.
    %
    %      (-1)   (0)   (1)
    %         1----3----2
    %

    properties (SetAccess = protected, GetAccess = public)    
        n_sides = 2;                    % no. "sides" (points in 1D)
        side_dof = [1; 2];              % local dofs of the "sides"
        side_type = 'Point';            % sides are points
        quad = ...                      % instance of Gauss quadrature class
            mFEM.Gauss('Order', 2, 'Type', 'line');
    end
    
    methods     
        function obj = Line3(id, nodes, varargin)
           % LINE3 Class constructor; calls base class constructor
           
           % Special case
           if ~all(size(nodes) == [2,1]) && ~all(size(nodes) == [3,1]);
                error('Line3:Line3','Nodes not specified correctly; expected a [2x1] or [3x1] array, but recieved a [%dx%d] array.', size(nodes,1), size(nodes,2));
           end
           
           % Test that nodes is sized correctly
           if all(size(nodes) == [2,1]);
               nodes(3,:)= mean(nodes,1);
           end
           
           % Call the base class constructor
           obj = obj@mFEM.Element(id, nodes, varargin{:}); 
           
           % Set the node plot order
           obj.node_plot_order = [1,3,2];
        end
        
        % Define the size function
        function L = size(obj)
        	L = norm(diff(obj.nodes,1));
        end
    end
    
    methods (Access = protected)          
        % Define the inherited abstract methods (protected)

        function N = basis(~, xi)
            % Returns a row vector of local shape functions
            N(1) = 1/2*xi*(xi-1);
            N(3) = 1-xi^2;
            N(2) = 1/2*xi*(1 + xi);
        end

        function B = grad_basis(obj, xi) 
            % Gradient of shape functions
            B = inv(obj.jacobian(xi)) * obj.local_grad_basis(xi);
        end
             
        function J = jacobian(obj, xi)
            % Returns the jacobian matrix
            J = obj.local_grad_basis(xi)*obj.nodes;                 
        end
        
        function GN = local_grad_basis(~, xi)
            GN = [xi-1/2, xi+1/2, -2*xi];
        end
    end
end