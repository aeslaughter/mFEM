classdef Quad4 < mFEM.Element
    %Quad4 4-node quadrilateral element
    %
    %   (-1,1)    (3)    (1,1)
    %         4---------3
    %         |         |
    %      (4)|         |(2)
    %         |         |
    %         1---------2
    %  (-1,-1)    (1)    (1,-1)
    %

    % Define the inherited abstract properties
    properties (SetAccess = protected, GetAccess = public)
        n_sides = 4;                        % no. of sides
        side_dof = [1,2; 2,3; 3,4; 4,1];    % define the side dofs 
        side_type = 'Truss2';               % side is 2-node truss element
        quad = mFEM.Gauss(2, 'quad');       % quadrature rules
     end
    
    % Define the Quad4 constructor
    methods
        function obj = Quad4(id, nodes, varargin)
           % Class constructor; calls base class constructor
           
           % Test that nodes is sized correctly
           if ~all(size(nodes) == [4,2]);
                error('Quad4:Quad4','Nodes not specified correctly; expected a [4x2] array, but recieved a [%d x %d] array.', size(nodes,1), size(nodes,2));
           end
           
           % Call the base class constructor
           obj = obj@mFEM.Element(id, nodes, varargin{:}); 
        end
    end
    
    % Define the inherited abstract methods (protected)
    methods (Access = protected)      
        function N = basis(~, xi, eta)
            % Returns a row vector of local shape functions
            N(1) = 1/4*(1-xi)*(1-eta);
            N(2) = 1/4*(1+xi)*(1-eta);
            N(3) = 1/4*(1+xi)*(1+eta);
            N(4) = 1/4*(1-xi)*(1+eta);
        end
        
        function GN = local_grad_basis(~, xi, eta)
            % Returns gradient, in xi and eta, of the shape functions
            GN = 1/4*[eta-1, 1-eta, 1+eta, -eta-1;
                      xi-1, -xi-1, 1+xi, 1-xi];
        end
        
        function B = grad_basis(obj, xi, eta) 
            % Gradient of shape functions
            B = inv(obj.jacobian(xi, eta)) * obj.local_grad_basis(xi, eta);
        end
        
        function J = jacobian(obj, xi, eta)
            % Returns the jacobian matrix  
            J = obj.local_grad_basis(xi, eta)*obj.nodes;                 
        end
    end
end