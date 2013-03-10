classdef Line2 < mFEM.elements.base.Element
    %LINE2 A 2-node, 1D linear element.
    %
    %      (-1)   (1)   (1)
    %         1---------2
    %
    %----------------------------------------------------------------------
    %  mFEM: An Object-Oriented MATLAB Finite Element Library
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
    %  along with this program.  If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
    
    properties (Constant)
        side_ids = [1; 2];              % local dofs of "sides"
        side_type = 'Point';            % sides are Point
        quad = mFEM.Gauss(1, 'line');   % Gauss quadrature class
               
    end
    
    methods     
        function obj = Line2(varargin)
            %LINE2 Class constructor; calls base class constructor
            %
            % Syntax
            %   Line2(id, nodes)
            %
            % Description
            %   see mFEM.Element
            %
            % See Also mFEM.Element
            obj = obj@mFEM.elements.base.Element(varargin{:}); 
        end
        
        % Define the size function
%         function L = size(obj)
%         	L = norm(diff(obj.nodes,1));
%         end
    end
    
    methods (Access = protected)          
%         function N = basis(~, xi)
%             % Returns a row vector of local shape functions
%             N(1) = 1/2*(1 - xi);
%             N(2) = 1/2*(1 + xi);
%         end
% 
%         function B = gradBasis(obj, varargin) 
%             % Gradient of shape functions
%             B = inv(obj.jacobian()) * obj.localGradBasis();
%         end
%              
%         function J = jacobian(obj, varargin)
%             % Returns the jacobian matrix (1/2 the length)
%             J = 1/2 * obj.size();                 
%         end
%         
%         function GN = localGradBasis(obj, varargin)
%             % Returns shape function derivatives in terms of xi
%             GN = [-1/2, 1/2];
%         end
    end
    
    methods (Static)
        function node_map = buildNodeMap(x0,x1,xn)
            spmd
                node_map = codistributed((x0:(x1-x0)/xn:x1)');
            end
        end

        function elem_map = buildElementMap(~,~,~,xn)
            spmd
                id = 1:xn+1;
                elem_map = zeros(xn,2,'uint32');
                for i = 1:xn;
                    elem_map(i,:) = [id(i),id(i+1)];
                end
                codist = codistributor1d(1,codistributor1d.unsetPartition,size(elem_map));
                elem_map = codistributed(elem_map,codist); 
            end
        end
    end
    
end