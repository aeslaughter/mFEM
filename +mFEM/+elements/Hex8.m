classdef Hex8 < mFEM.base.Element
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
    properties (Constant)
        side_ids =...                      % define the side dofs 
            [1,2,3,4; 5,6,7,8; 1,2,6,5; 2,3,7,6; 3,4,8,7; 4,1,5,8];
        side_type = 'Quad4';                % side is 4-node Quad element
        quad = mFEM.Gauss(2,'hex');         % quadrature rules
        n_dim = 3;
        n_nodes = 8;
        node_plot_order = [1,2,3,4,8,7,6,5];
     end
    
    % Define the Quad4 constructor
    methods
        function obj = Hex8(varargin)
           %HEX8 Class constructor; calls base class constructor
           
           % Call the base class constructor
           obj = obj@mFEM.base.Element(varargin{:}); 
           
           % Set the node plotting order
%            obj.node_plot_order = [1,2,3,4,8,7,6,5];  
        end
    end
    
    % Define the inherited abstract methods (protected)
    methods (Access = protected)      
        function N = basis(~, x)
            %BASIS returns the element basis functions
            
            % Convert x to xi,eta,zeta for convience
            xi = x(1); eta = x(2); zeta = x(3);
            
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
        
        function GN = localGradBasis(~, x)
            %LOCALGRADBASIS returns the basis functions derivatives
            
            % Convert x to xi,eta,zeta for convience
            xi = x(1); eta = x(2); zeta = x(3);
            
            % Define the derivatives (see Hex8_deriv.m)
            GN(1,:) = [ -((eta - 1)*(zeta - 1))/8, ((eta - 1)*(zeta - 1))/8, -((eta + 1)*(zeta - 1))/8, ((eta + 1)*(zeta - 1))/8, ((eta - 1)*(zeta + 1))/8, -((eta - 1)*(zeta + 1))/8, ((eta + 1)*(zeta + 1))/8, -((eta + 1)*(zeta + 1))/8];
            GN(2,:) = [ -(xi/8 - 1/8)*(zeta - 1), (xi/8 + 1/8)*(zeta - 1), -(xi/8 + 1/8)*(zeta - 1), (xi/8 - 1/8)*(zeta - 1), (xi/8 - 1/8)*(zeta + 1), -(xi/8 + 1/8)*(zeta + 1), (xi/8 + 1/8)*(zeta + 1), -(xi/8 - 1/8)*(zeta + 1)];
            GN(3,:) = [ -(xi/8 - 1/8)*(eta - 1), (xi/8 + 1/8)*(eta - 1), -(xi/8 + 1/8)*(eta + 1), (xi/8 - 1/8)*(eta + 1), (xi/8 - 1/8)*(eta - 1), -(xi/8 + 1/8)*(eta - 1), (xi/8 + 1/8)*(eta + 1), -(xi/8 - 1/8)*(eta + 1)];
        end
        
        function B = gradBasis(obj, x) 
            %GRADBASIS Gradient of shape functions in x,y,z
            B = inv(obj.jacobian(x)) * obj.localGradBasis(x);
        end
        
        function J = jacobian(obj, x)
            %JACOBIAN Returns the jacobian matrix  
            J = obj.localGradBasis(x)*obj.nodes;               
        end
    end
    
    methods (Static)
        function node_map = buildNodeMap(x0,x1,y0,y1,z0,z1,xn,yn,zn)
            spmd
                x = codistributed(x0:(x1-x0)/xn:x1);
                y = codistributed(y0:(y1-y0)/yn:y1);
                z = codistributed(z0:(z1-z0)/zn:z1);
                [X,Y,Z] = ndgrid(x,y,z); 
                node_map = [reshape(X,numel(X),1),...
                            reshape(Y,numel(Y),1),...
                            reshape(Z,numel(Z),1)];
            
                if numlabs == 1;
                    codist = codistributor1d(1,codistributor1d.unsetPartition,size(node_map));
                    node_map = redistribute(node_map,codist);
                end
            end
        end
            
        function elem_map = buildElementMap(~,~,~,~,~,~,xn,yn,zn)
            spmd
                n = (xn+1)*(yn+1)*(zn+1);
                id = reshape(1:n,xn+1,yn+1,zn+1);
                l = 0;
                elem_map = zeros(xn*yn*zn,8,'uint32');
                for k = 1:zn;
                    for j = 1:yn;
                        for i = 1:xn;
                            l = l + 1;
                            elem_map(l,:) = [id(i,j,k), id(i+1,j,k), id(i+1,j,k+1), id(i,j,k+1),...
                                             id(i,j+1,k), id(i+1,j+1,k), id(i+1,j+1,k+1), id(i,j+1,k+1)];
                        end
                    end
                end
                
                codist = codistributor1d(1,codistributor1d.unsetPartition,size(elem_map));
                elem_map = codistributed(elem_map,codist);
            end
        end
    end  
end