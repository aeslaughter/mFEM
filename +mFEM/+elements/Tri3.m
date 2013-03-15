classdef Tri3 < mFEM.base.Element
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
    
    % Define the inherited abstract properties
    properties (Constant)
        side_ids = [1,2; 2,3; 3,1]; % define the side dofs 
        side_type = 'Line2';        % 2-node Truss element for side
        quad = mFEM.Gauss(3, 'tri');% 3-point triangular quadrature 
        n_dim = 2;                  % 2D space
        n_nodes = 3;                % no. of nodes
    end
    
    % Define the Quad4 constructor
    methods
        function obj = Tri3(varargin)
           % Class constructor; calls base class constructor
           obj = obj@mFEM.base.Element(varargin{:});
        end
    end
    
    % Define the inherited abstract methods (protected)
    methods (Access = protected)    
        
        function J = jacobian(obj, ~)
            %JACOBIAN Define the Jacobain (see Fish p. 180)
            %  Note, the Jacobian is a constant.
           
            % Extract node coordinates
            nodes = obj.nodes.getCoord();
              
            % Define short-hand for difference between two points
            x = @(i,j) nodes(i,1) - nodes(j,1);
            y = @(i,j) nodes(i,2) - nodes(j,2);
            
            % Define the Jacobian matrix
            J = [x(1,3), y(1,3); x(2,3), y(2,3)];
        end
        
        function N = basis(~, x)
            %BASIS Returns a row vector of local shape functions
            
            % Shape function vector
            N(1) = x(1);
            N(2) = x(2);
            N(3) = 1 - sum(x);
          end

        function B = gradBasis(obj, ~) 
            %GRADBASIS Gradient of shape functions (Fish, p. 174)
            %   Note, the gradient of N is constant.
           
            % Extract node coordinates
            nodes = obj.nodes.getCoord();
            
            % Define short-hand for difference between two points
            x = @(i,j) nodes(i,1) - nodes(j,1);
            y = @(i,j) nodes(i,2) - nodes(j,2);
        
            % Define twice the area of element
            M = [ones(obj.n_nodes,1), nodes]; % 2A

            % The shape function derivatives
            B = 1/det(M)*[y(2,3), y(3,1), y(1,2);
                          x(3,2), x(1,3), x(2,1)];
        end
        
        function GN = localGradBasis(obj, ~)
            error('Tri3:local_grad_basis', 'Function not defined for the %s element, the B matrix is computed directly.', class(obj));
        end
    end
end