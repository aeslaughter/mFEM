classdef Point < mFEM.elements.base.Element
    %POINT A 1-node "element", for use as side element of 1D elements.
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
    
    properties (Constant, GetAccess = public)    
        n_sides = [];                   % no. "sides" (nodes are sides in 1D)
        side_dof = [];                  % local dofs of the "sides"
        side_type = '';                 % 1D elements do not have side elements
        quad = [];                      % Instance of Gauss quadrature class
    end
    
    methods     
        function obj = Point(id, nodes, varargin)
           % POINT Class constructor; calls base class constructor

           % Call the base class constructor
           obj = obj@mFEM.elements.base.Element(id, nodes, varargin{:}); 
        end
    end
    
    methods (Access = protected)          
        function N = basis(obj,varargin)
            % Returns a row vector of local shape functions
            N = 1;
        end

        function B = gradBasis(obj, varargin) 
            % Gradient of shape functions
            error('Point:gradBasis', 'Not defined for Point element');
            B = inv(obj.jacobian()) * obj.local_grad_basis;
        end
             
        function J = jacobian(obj, varargin)
            % Returns the jacobian matrix  
            error('Point:jacobian', 'Not defined for Point element'); 
        end
        
        function GN = localGradBasis(obj, varargin)
            error('Point:localGradBasis', 'Not defined for Point element');
        end
    end
end