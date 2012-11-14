classdef Truss < mFEM.Element
    % A 2-node Truss element, it may be located in 1D, 2D, or 3D space.
    %
    %      (-1)   (1)   (1)
    %         1---------2
    %
    % This is a special case...
    %

    % Define the inherited abstract properties
    properties (SetAccess = protected, GetAccess = public)
        n_sides = 2;                 % no. "sides" (nodes are sides in 1D)
        side_dof = [1; 2];           % local dofs of the "sides"
        side_type = 'Point';         % 1D elements have points on sides
        quad = []                    % quadrature not necessary
    end
    
    methods  
        % Define the Truss constructor
        function obj = Truss(id, nodes, varargin)
           % Class constructor; calls base class constructor
           
           % Test that nodes is sized correctly
           if ~all(size(nodes) == [2,2]) && ~all(size(nodes) == [2,3]) ;
               error('Truss:Truss','Nodes not specified correctly; expected a [2x2] or [2x3] array, but recieved a [%dx%d] array.', size(nodes,1), size(nodes,2));
           end
           
           % Call the base class constructor
           obj = obj@mFEM.Element(id, nodes, 'Space', 2); 
           
           % Set the local dimensionality, only needed if the number of 
           % element local coordinates (xi,eta,...) are different from the 
           % number of spatial coordinates (x,y,...).
           %obj.local_n_dim = 1;
           
           % Set the number of dofs per node on the element
           %obj.n_dof_node = obj.n_dim; 
        end
        
        % Define the size function
        function L = size(obj)
        	L = norm(diff(obj.nodes,1));
        end
        
        function K = stiffness(obj, varargin)
           
            N = [1,0,-1,0];
                       
            d = diff(obj.nodes)/obj.size();
            c = d(1); s = d(2);

            T = zeros(4,4);
            T(1:2,1:2) = [c,s;-s,c];
            T(3:4,3:4) = [c,s;-s,c];
            
            
            K = T'*(N'*N)*T;
        end
        
%         function T = transformation(obj, varargin)
%            %TRANSFORMATION Outputs the transformation matrix, T
% 
%         end
%         function N = shape(obj, varargin)
%            N = zeros(1,obj.n_nodes*obj.n_dim);
%            N(1:obj.n_dim:end) = [1,-1];     
%         end
    end
    
    methods (Access = protected)  
        function basis(obj, varargin)
            % Returns a row vector of local shape functions
            error('Truss:basis','The gradient of the basis functions is not defined for the Truss element.');   
        end

        function grad_basis(~, varargin) 
            % Gradient of shape functions
            error('Truss:grad_basis','The gradient of the basis functions is not defined for the Truss element.');   
        end
             
        function jacobian(~, varargin)
            % Returns the jacobian matrix (1/2 the length)
            error('Truss:jacobian','The jacabian matrix is not defined for the Truss element.');   
        end
        
        function local_grad_basis(~, varargin)
            % Does nothing for the Truss dlement
            error('Truss:local_grad_basis','The gradient in local coordinte system is not defined for the Truss element.');   
        end
    end
end