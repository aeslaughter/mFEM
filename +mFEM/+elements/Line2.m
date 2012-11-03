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
            %   Line2(id, nodes) creates an element given, where id is a
            %   unique identification number for this element and nodes is 
            %   a matrix of node coordinates (global) that should be 
            %   arranged as column matrix (no. nodes x no. dims).
            %
            %   Line2(id, nodes, 'PropertyName', PropertyValue, ...) 
            %   allows the user to customize the behavior of the element, 
            %   the available properties are listed below.
            %
            % Element Property Descriptions
            %   space
            %       'scalar' | 'vector' | integer
            %       Allows the type of FEM space to be set: scalar sets the 
            %       number of dofs per node to 1, vector  sets it to the 
            %       no. of space dimension, and  specifing a number sets it
            %       to that value.
            %
            % See Also mFEM.Element
           
            % Test that nodes is sized correctly
            if ~all(size(nodes) == [2,1]);
                error('Line2:Line2','Nodes not specified correctly; expected a [2x1] array, but recieved a [%dx%d] array.', size(nodes,1), size(nodes,2));
            end

            % Call the base class constructor
            obj = obj@mFEM.Element(id, nodes, varargin{:}); 
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
            % Returns the jacobian matrix  
            J = obj.local_grad_basis()*obj.nodes;                 
        end
        
        function GN = local_grad_basis(obj, varargin)
            % Returns shape function derivatives in terms of xi
            GN = [-1/2, 1/2];
        end
    end
end