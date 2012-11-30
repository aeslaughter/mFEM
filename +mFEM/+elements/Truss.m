classdef Truss < mFEM.base.ElementCore
    % A 2-node Truss element, it may be located in 1D, 2D, or 3D space.
    %
    %      (-1)   (1)   (1)
    %         1---------2
    %
    % This is a special case...
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
           obj = obj@mFEM.base.ElementCore(id, nodes, 'Space', 2); 
           
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
        

    end
    
  
end