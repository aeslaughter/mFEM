classdef Line3 < mFEM.Element
    %LINE3 A 3-node, 1D quadradtic element.
    %
    %      (-1)   (0)   (1)
    %         1----2----3
    %

    properties (SetAccess = protected, GetAccess = public)    
        n_sides = 2;                    % no. "sides" (points in 1D)
        side_dof = [1; 3];              % local dofs of the "sides"
        side_type = 'Point';            % Sides are points
        quad = mFEM.Gauss(2,'line');    % Instance of Gauss quadrature class
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
               mid = mean(nodes,1);
               nodes(3,:) = nodes(2,:);
               nodes(2,:) = mid;
           end
           
           % Call the base class constructor
           obj = obj@mFEM.Element(id, nodes, varargin{:}); 
        end
    end
    
    methods (Access = protected)          
        % Define the inherited abstract methods (protected)

        function N = basis(~, xi)
            % Returns a row vector of local shape functions
            
            %%%% ADD SHAPE FUNCTIONS HERE %%%%
            
        end

        function B = grad_basis(obj, xi) 
            % Gradient of shape functions
            
            %%%% ADD GRADIENT SHAPE FUNCTIONS (dN/dx) HERE %%%%
            
        end
             
        function J = jacobian(obj, xi)
            % Returns the jacobian matrix
            
            %%%% ADD CALCULATION OF JACOBIAN HERE %%%%
               
        end
        
        function GN = local_grad_basis(~, xi)
            % Return shape function derivatives, in terms of xi
            
            %%%% ADD SHAPE FUNCTION DERIVATIVES (dN/dxi) HERE %%%%

        end
    end
end