classdef Quad4 < mFEM.elements.base.Element
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
        n_nodes = uint32(4);
        n_sides = uint32(4);
        side_ids = [1,2; 2,3; 3,4; 4,1];    % define the side dofs 
%         side_type = 'Line2';                % side is 2-node line element
%         quad = ...                          % quadrature rules
%             mFEM.Gauss('order', 2, 'type', 'quad');       
     end
    
    % Define the Quad4 constructor
    methods
        function obj = Quad4(varargin)
           % Class constructor; calls base class constructor
           obj = obj@mFEM.elements.base.Element(varargin{:}); 
        end
    end
    
    % Define the inherited abstract methods (protected)
    methods (Access = protected)      
%         function N = basis(~, x)
%             %BASIS Returns a row vector of local shape functions
% 
%             % Define xi and eta from vector input
%             xi = x(1);
%             eta = x(2);
%             
%             % Compute the shape function vector
%             N(1) = 1/4*(1-xi)*(1-eta);
%             N(2) = 1/4*(1+xi)*(1-eta);
%             N(3) = 1/4*(1+xi)*(1+eta);
%             N(4) = 1/4*(1-xi)*(1+eta);
%         end
%         
%         function GN = localGradBasis(~, x)
%             %LOCALGRADBASIS gradient, in xi and eta, of shape functions
%             
%             % Define xi and eta from vector input
%             xi = x(1);
%             eta = x(2);
%             
%             % Compute gradient
%             GN = 1/4*[eta-1, 1-eta, 1+eta, -eta-1;
%                       xi-1, -xi-1, 1+xi, 1-xi];
%         end
%         
%         function B = gradBasis(obj, x) 
%             %GRADBASIS Gradient of shape functions
%                         
%             % Compute the gradient of bais in x,y
%             B = inv(obj.jacobian(x)) * obj.localGradBasis(x);
%         end
%         
%         function J = jacobian(obj, x)
%             %JACOBIAN Returns the jacobian matrix  
%                         
%             % Compute the Jacobian
%             J = obj.localGradBasis(x)*obj.nodes;                 
%         end
    end
    
    methods (Static)%(Static, Access = ?mFEM.Mesh)
        function p = plot(vert, face)
            p = patch('Vertices',vert,'Faces',face);
        end
        
        function grid(x0,x1,y0,y1,xn,yn)
            
            tic;
            [nodes, node_map] = buildGrid(x0,x1,y0,y1,xn,yn);
            toc;

%             n_elem = xn*yn;
%             n_nodes = length(nodes);
%             elem_id = reshape(1:n_nodes,xn+1,yn+1);
% 
%             elem_map = zeros(n_elem,4,'uint32');
%             id = 0;
%             for j = 1:yn;
%                 for i = 1:xn;
%                     id = id + 1;
%                     idx = [elem_id(i,j), elem_id(i+1,j),...
%                            elem_id(i+1,j+1), elem_id(i,j+1)];
%                     elem_map(id,:) = idx;
%                 end
%             end
%             
%             nodes = gather(nodes);
%             elements = cell(n_elem,1);
%             for i = 1:n_elem;
%                 no = nodes(elem_map(i,:));
%                 elements{i} = mFEM.elements.Quad4(i,no);
%             end
%             
%             elements = distributed(elements);
            
        end
    end    
end