classdef Truss < mFEM.base.ElementDirect
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
           obj = obj@mFEM.base.ElementDirect(id, nodes, 'Space', 2); 
           
        end
        
        % Define the size function
        function L = size(obj)
        	L = norm(diff(obj.nodes,1));
        end
        
        function Ke = stiffness(obj, varargin)
           
            N = [1,0,-1,0];
                       
            d = diff(obj.nodes)/obj.size();
            c = d(1); s = d(2);

            T = zeros(4,4);
            T(1:2,1:2) = [c,s;-s,c];
            T(3:4,3:4) = [c,s;-s,c];
            
            
            Ke = T'*(N'*N)*T;
        end
        
        function fe = force(obj, varargin)
            
        end

    end
    
  
end