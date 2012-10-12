classdef Tri3 < mFEM.Element
    %Tri3 3-node triangle element
    %
    % This element is a Quad4 with nodes 3 and 4 shared.
    %
    %   (-1,1)   
    %         (3)
    %         |\       
    %      (3)| \(2)      
    %         |  \       
    %         1---2
    %  (-1,-1) (1)  (1,-1)
    %

    % Define the inherited abstract properties
    properties (SetAccess = protected, GetAccess = public)
        n_shape = 4; % no. of shape functions
        n_sides = 3; % no. of sides
        side_dof = [1,2; 2,3; 3,1];         % define the side dofs 
        side_defn = [2,-1; 1,1; 1,-1];      % xi,eta definitions for sides
        collapsed_nodes = [4,3];            % node 4 is collapsed to node 3
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
        function N = basis(~, xi, eta)
            % Returns a row vector of local shape functions
            N(1) = 1/4*(1-xi)*(1-eta);
            N(2) = 1/4*(1+xi)*(1-eta);
            N(3) = 1/4*(1+xi)*(1+eta);
        %    N(4) = 1/4*(1-xi)*(1+eta);
        end

        function GN = grad_basis(~, xi, eta) 
            % Gradient of shape functions
            GN = 1/4*[eta-1, 1-eta, 1+eta, -eta-1;
                      xi-1, -xi-1, 1+xi, 1-xi];
        end
    end
end