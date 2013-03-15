classdef Truss < mFEM.base.Element
    %TRUSS A 2-node Truss element, it may be in 1D, 2D, or 3D space.
    %
    %      (-1)   (1)   (1)
    %         1---------2
    %
    % This is a special case...
    %
    %----------------------------------------------------------------------
    %  mFEM: A Parallel, Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2013 Andrew E Slaughter
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
    %  along with this program. If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------

    % Define the inherited abstract properties
    properties (Constant)
        side_ids = [1; 2];           % local dofs of the "sides"
        side_type = 'Point';         % 1D elements have points on sides
        quad = []                    % quadrature not necessary
        n_dim = 2;                   % no. of space dimensions.
        n_nodes = 2;                 % no. of nodes
    end
    
    methods  
        % Define the Truss constructor
        function obj = Truss(varargin)
           % Class constructor; calls base class constructor
           obj = obj@mFEM.base.Element(varargin{:}); 
        end
        
        % Define the size function
        function L = size(obj)
        	L = norm(diff(obj.nodes.getCoord(),1));
        end

        function Ke = stiffness(obj, varargin)
           
            N = [1,0,-1,0];
                       
            d = diff(obj.nodes.getCoord())/obj.size();
            c = d(1); s = d(2);

            T = zeros(4,4);
            T(1:2,1:2) = [c,s;-s,c];
            T(3:4,3:4) = [c,s;-s,c];

            Ke = T'*(N'*N)*T;
        end
    end
            
    methods (Access = protected)       
        function basis(varargin)
            error('Truss:basis','Not implemented for the Truss element.');
        end
        
        function localGradBasis(varargin)
            error('Truss:localGradBasis','Not implemented for the Truss element.');
        end
        
        function gradBasis(varargin) 
            error('Truss:gradBasis','Not implemented for the Truss element.');
        end
        
        function jacobian(varargin)
            error('Truss:jacobian','Not implemented for the Truss element.');
        end
    end
end