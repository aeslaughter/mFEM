classdef Quad4 < Element
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
        n_shape = 4; % no. of shape functions
        n_sides = 4; % no. of sides
        side_dof = [1,2; 2,3; 3,4; 4,1];  
     end
    
    % Define the Quad4 constructor
    methods
        function obj = Quad4(id, nodes, varargin)
           % Class constructor; calls base class constructor
           obj = obj@Element(id, nodes, varargin{:}); 
        end
        
        % (...automate this...)
        function side = build_side(obj, id)
            % Creates a Linear2 element for the specified side

            % Extract the nodes for the side
            idx = obj.side_dof(id, :);
            nodes = obj.nodes(idx,:);       

            % Build a local coordinate system (0 to length of side)
            n = [0; norm(nodes(1,:) - nodes(2,:))]; 

            % Create the element
            side = Linear2(id, n, obj.space);
            
            % Add boundary_ids to the created side
            side.boundary_id = obj.side(id).boundary_id;
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

        function GN = grad_basis(~, xi, eta) 
            % Gradient of shape functions
            GN = 1/4*[eta-1, 1-eta, 1+eta, -eta-1;
                      xi-1, -xi-1, 1+xi, 1-xi];
        end
        

    end
end