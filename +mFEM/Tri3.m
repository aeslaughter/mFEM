classdef Tri3 < mFEM.Element
    %Tri3 3-node triangle element
    % 
    %        xi2  
    %         ^
    %         |  
    %
    %         2
    %         |\            On side 1: xi1 + xi2 = 1; 
    %      (2)| \(1)        On node 3: xi3 = 1;
    %         |  \       
    %         3---1  ---> xi1
    %          (3) 

    % Define the inherited abstract properties
    properties (SetAccess = protected, GetAccess = public)
        n_shape = 3; % no. of shape functoins
        n_sides = 3; % no. of sides
        side_dof = [1,2; 2,3; 3,1]; % define the side dofs 
        side_type = 'Linear2';
        side_defn = [2,0; 1, 0; NaN, NaN];      % xi1,xi2 definitions for sides
     end
    
    % Define the Quad4 constructor
    methods
        function obj = Tri3(id, nodes, varargin)
           % Class constructor; calls base class constructor
           
           % Test that nodes is sized correctly
           if ~all(size(nodes) == [3,2]);
                error('Tri3:Tri3','Nodes not specified correctly; expected a [3 x 2] array, but recieved a [%d x %d] array.', size(nodes,1), size(nodes,2));
           end
           
           % Call the base class constructor
           obj = obj@mFEM.Element(id, nodes, varargin{:});
        end
    end
    
    % Define the inherited abstract methods (protected)
    methods (Access = protected)    
        
        function J = jacobian(obj, varargin)
            % Define the Jacobain (see Fish p. 180)
            
            % The Jacobian is constant, produce a warning for agruements
            if nargin > 1;
                warning('Tri3:jacobian', 'The Jacobian for the Tri3 element is constant, thus no spatial coordinates are needed.');
            end
            
            % Define short-hand for difference between two points
            x = @(i,j) obj.nodes(i,1) - obj.nodes(j,1);
            y = @(i,j) obj.nodes(i,2) - obj.nodes(j,2);
            
            % Define the Jacobian matrix
            J = [x(1,3), y(1,3); x(2,3), y(2,3)];
        end
        
        function N = basis(~, xi1, xi2)
            % Returns a row vector of local shape functions
            N(1) = xi1;
            N(2) = xi2;
            N(3) = 1 - xi2 - xi1;
          end

        function B = grad_basis(obj, varargin) 
            % Gradient of shape functions (Fish, p. 174)
            
            % The Jacobian is constant, produce a warning for agruements
            if nargin > 1;
                warning('Tri3:grad_basis', 'The shape functin dervatives for the Tri3 element are constant, thus no spatial coordinates are needed.');
            end
            % Define short-hand for difference between two points
            x = @(i,j) obj.nodes(i,1) - obj.nodes(j,1);
            y = @(i,j) obj.nodes(i,2) - obj.nodes(j,2);
        
            % Define twice the area of element
            M = [ones(obj.n_nodes,1), obj.nodes]; % 2A

            % The shape function derivatives
            B = 1/det(M)*[y(2,3), y(3,1), y(1,2);
                           x(3,2), x(1,3), x(2,1)];
        end
    end
end