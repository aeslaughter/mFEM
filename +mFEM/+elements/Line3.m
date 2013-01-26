classdef Line3 < mFEM.elements.base.Element
    %LINE3 A 3-node, 1D quadradtic element.
    %
    %      (-1)   (0)   (1)
    %         1----3----2
    %----------------------------------------------------------------------
    %  mFEM: An Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2012 Andrew E Slaughter
    % 
    %  This program is free software: you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation, either version 3 of the License, or
    %  (at your option) any later version.
    % 
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    % 
    %  You should have received a copy of the GNU General Public License
    %  along with this program.  If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
    
    properties (SetAccess = protected, GetAccess = public)    
        n_sides = 2;                    % no. "sides" (points in 1D)
        side_dof = [1; 2];              % local dofs of the "sides"
        side_type = 'Point';            % sides are points
        quad = ...                      % instance of Gauss quadrature class
            mFEM.Gauss('Order', 2, 'Type', 'line');
    end
    
    methods     
        function obj = Line3(id, nodes, varargin)
           % LINE3 Class constructor; calls base class constructor
           
           % Special case
           if ~all(size(nodes) == [2,1]) && ~all(size(nodes) == [3,1]) && ~all(size(nodes) == [3,2]);;
                error('Line3:Line3','Nodes not specified correctly; expected a [2x1], [3x1], or [3x2] array, but recieved a [%dx%d] array.', size(nodes,1), size(nodes,2));
           end
           
           % Test that nodes is sized correctly
           if all(size(nodes) == [2,1]);
               nodes(3,:)= mean(nodes,1);
           end
           
           % Call the base class constructor
           obj = obj@mFEM.elements.base.Element(id, nodes, varargin{:}); 
           
           % Set the node plot order
           obj.node_plot_order = [1,3,2];
        end
        
        % Define the size function
        function L = size(obj)
        	L = norm(diff(obj.nodes,1));
        end
    end
    
    methods (Access = protected)          
        % Define the inherited abstract methods (protected)

        function N = basis(~, xi)
            % Returns a row vector of local shape functions
            N(1) = 1/2*xi*(xi-1);
            N(3) = 1-xi^2;
            N(2) = 1/2*xi*(1 + xi);
        end

        function B = gradBasis(obj, xi) 
            % Gradient of shape functions
            B = inv(obj.jacobian(xi)) * obj.localGradBasis(xi);
        end
             
        function J = jacobian(obj, xi)
            % Returns the jacobian matrix
            J = obj.localGradBasis(xi)*obj.nodes;                 
        end
        
        function GN = localGradBasis(~, xi)
            GN = [xi-1/2, xi+1/2, -2*xi];
        end
    end
end