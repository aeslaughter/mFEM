classdef Point < mFEM.base.Element
    %POINT A 1-node "element", for use as side element of 1D elements.

    properties (SetAccess = protected, GetAccess = public)    
        n_sides = [];                   % no. "sides" (nodes are sides in 1D)
        side_dof = [];                  % local dofs of the "sides"
        side_type = '';                 % 1D elements do not have side elements
        quad = [];                      % Instance of Gauss quadrature class
    end
    
    methods     
        function obj = Point(id, nodes, varargin)
           % POINT Class constructor; calls base class constructor

           % Call the base class constructor
           obj = obj@mFEM.base.Element(id, nodes, varargin{:}); 
        end
    end
    
    methods (Access = protected)          
        function N = basis(~)
            % Returns a row vector of local shape functions
            N = 1;
        end

        function B = grad_basis(obj, varargin) 
            % Gradient of shape functions
            error('Point:grad_basis', 'Not defined for Point element');
            B = inv(obj.jacobian()) * obj.local_grad_basis;
        end
             
        function J = jacobian(obj, varargin)
            % Returns the jacobian matrix  
            error('Point:jacobian', 'Not defined for Point element'); 
        end
        
        function GN = local_grad_basis(obj, varargin)
            error('Point:local_grad_basis', 'Not defined for Point element');
        end
    end
end