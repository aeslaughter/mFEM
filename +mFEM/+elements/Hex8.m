classdef Hex8 < mFEM.elements.base.Element
    %Hex8 8-node hexahedron element
    %
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
    properties (SetAccess = protected, GetAccess = public)
        n_sides = 6;                       % no. of sides
        side_dof =...                     % define the side dofs 
            [1,2,3,4; 5,6,7,8; 1,2,6,5; 2,3,7,6; 3,4,8,7; 4,1,5,8];
        side_type = 'Quad4';                % side is 4-node Quad element
        quad = ...                          % quadrature rules
            mFEM.Gauss('order', 2, 'type', 'hex');       
     end
    
    % Define the Quad4 constructor
    methods
        function obj = Hex8(id, nodes, varargin)
           % Class constructor; calls base class constructor
           
           % Test that nodes is sized correctly
           if ~all(size(nodes) == [8,3]);
                error('Hex8:Hex8','Nodes not specified correctly; expected a [8x3] array, but recieved a [%d x %d] array.', size(nodes,1), size(nodes,2));
           end
           
           % Call the base class constructor
           obj = obj@mFEM.elements.base.Element(id, nodes, varargin{:}); 
           
           % Set the node plotting order
           obj.node_plot_order = [1,2,3,4,8,7,6,5]; 
           
        end
    end
    
    % Define the inherited abstract methods (protected)
    methods (Access = protected)      
        function N = basis(~, xi, eta, zeta)
            % Returns a row vector of local shape functions
            N(1) = 1/8*(1-xi)*(1-eta)*(1-zeta);
            N(2) = 1/8*(1+xi)*(1-eta)*(1-zeta);
            N(3) = 1/8*(1+xi)*(1+eta)*(1-zeta);
            N(4) = 1/8*(1-xi)*(1+eta)*(1-zeta);
            N(5) = 1/8*(1-xi)*(1-eta)*(1+zeta);
            N(6) = 1/8*(1+xi)*(1-eta)*(1+zeta);
            N(7) = 1/8*(1+xi)*(1+eta)*(1+zeta);
            N(8) = 1/8*(1-xi)*(1+eta)*(1+zeta);
        end
        
        function GN = local_grad_basis(~, xi, eta, zeta)
            % Returns gradient, in xi and eta, of the shape functions
            GN(:,1) = 1/8*[(eta-1)*(zeta-1); (zeta-1)*(xi-1); (eta-1)*(xi-1)];
            GN(:,2) = 1/8*[(eta-1)*(zeta-1); (zeta-1)*(xi+1); (eta-1)*(xi+1)];
            GN(:,3) = 1/8*[(eta+1)*(zeta-1); (zeta-1)*(xi+1); (eta+1)*(xi+1)];
            GN(:,4) = 1/8*[(eta+1)*(zeta-1); (zeta-1)*(xi-1); (eta+1)*(xi-1)];
            GN(:,5) = 1/8*[(eta-1)*(zeta+1); (zeta+1)*(xi-1); (eta-1)*(xi-1)];
            GN(:,6) = 1/8*[(eta-1)*(zeta+1); (zeta+1)*(xi+1); (eta-1)*(xi+1)];
            GN(:,7) = 1/8*[(eta+1)*(zeta+1); (zeta+1)*(xi+1); (eta+1)*(xi+1)];
            GN(:,8) = 1/8*[(eta+1)*(zeta+1); (zeta+1)*(xi-1); (eta+1)*(xi-1)];
        end
        
        function B = grad_basis(obj, xi, eta,zeta) 
            % Gradient of shape functions
            B = inv(obj.jacobian(xi, eta, zeta)) * obj.local_grad_basis(xi, eta, zeta);
        end
        
        function J = jacobian(obj, xi, eta, zeta)
            % Returns the jacobian matrix  
            J = obj.local_grad_basis(xi, eta, zeta)*obj.nodes;                 
        end
    end
end