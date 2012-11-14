classdef Line2 < mFEM.Element
    %LINE2 A 2-node, 1D linear element.
    %
    %      (-1)   (1)   (1)
    %         1---------2
    %

    properties (SetAccess = protected, GetAccess = public)    
        n_sides = 2;                    % no. "sides" (points in 1D)
        side_dof = [1; 2];              % local dofs of the "sides"
        side_type = 'Point';            % sides are points
        quad = ...                      % instance of Gauss quadrature class
            mFEM.Gauss('Order', 1, 'Type', 'line');    
    end
    
    methods     
        function obj = Line2(id, nodes, varargin)
            %LINE2 Class constructor; calls base class constructor
            %
            % Syntax
            %   Line2(id, nodes)
            %   Line2(id, nodes, 'PropertyName', PropertyValue, ...)
            %
            % Description
            %   see mFEM.Element
            %
            % See Also mFEM.Element
           
            % Test that nodes is sized correctly
            if ~all(size(nodes) == [2,1]) && ~all(size(nodes) == [2,2]) && ~all(size(nodes) == [2,3]) ;
                error('Line2:Line2','Nodes not specified correctly; expected a [2x1], [2x2], or [2x3] array, but recieved a [%dx%d] array.', size(nodes,1), size(nodes,2));
            end

            % Call the base class constructor
            obj = obj@mFEM.Element(id, nodes, varargin{:}); 
        end
        
        % Define the size function
        function L = size(obj)
        	L = norm(diff(obj.nodes,1));
        end
    end
    
    methods (Access = protected)          
        function N = basis(~, xi)
            % Returns a row vector of local shape functions
            N(1) = 1/2*(1 - xi);
            N(2) = 1/2*(1 + xi);
        end

        function B = grad_basis(obj, varargin) 
            % Gradient of shape functions
            B = inv(obj.jacobian()) * obj.local_grad_basis();
        end
             
        function J = jacobian(obj, varargin)
            % Returns the jacobian matrix (1/2 the length)
            J = 1/2 * obj.size();                 
        end
        
        function GN = local_grad_basis(obj, varargin)
            % Returns shape function derivatives in terms of xi
            GN = [-1/2, 1/2];
        end
    end
end