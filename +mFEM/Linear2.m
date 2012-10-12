classdef Linear2 < mFEM.Element
    %Linear 2-node 1D linear element
    %
    %      (-1)   (1)   (1)
    %         1---------2
    %

    % Define the inherited abstract properties
    properties (SetAccess = protected, GetAccess = public)
        n_shape = 2;                % no. of shape functions
        n_sides = 2;                % no. "sides" (nodes are sides in 1D)
        side_dof = [1; 2];          % local dofs of the "sides"
        side_defn = [1,-1; 1,1];   % xi definitions for "sides"
    end
    
    % Define the Linear2 constructor
    methods 
        function obj = Linear2(id, nodes, varargin)
           % Class constructor; calls base class constructor
           obj = obj@mFEM.Element(id, nodes, varargin{:}); 
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