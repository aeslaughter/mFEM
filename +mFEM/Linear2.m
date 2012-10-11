classdef Linear2 < Element
    %Linear 2-node 1D linear element
    %
    %      (-1)   (1)   (1)
    %         1---------2
    %

    % Define the inherited abstract properties
    properties (SetAccess = protected, GetAccess = public)
        n_shape = 2; % no. of shape functions
        n_sides = [];
        side_elem = '';
        side_dof = [];
    end
    
    % Define the Quad4 constructor
    methods 
        function obj = Linear2(id, nodes, varargin)
           % Class constructor; calls base class constructor
           obj = obj@Element(id, nodes, varargin{:}); 
        end
        
        function build_side(~,~)
           % Class for building a side element
           error('ERROR: Side not defined for Linear2 element'); 
        end
    end
    
    % Define the inherited abstract methods (protected)
    methods (Access = protected)      
        function N = basis(~, xi)
            % Returns a row vector of local shape functions
            N(1) = 1/2*(1-xi);
            N(2) = 1/2*(1+xi);
        end

        function GN = grad_basis(~, ~) 
            % Gradient of shape functions
            GN = [-1/2, 1/2];
        end
    end
end