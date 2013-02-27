classdef Line2 < mFEM.cells.base.Cell
    %LINE2 A 2-node, 1D linear element.
    %
    %      (-1)   (1)   (1)
    %         1---------2
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
    
    properties (SetAccess = protected, GetAccess = public) 
        n_sides = 2;
        side_ids = [1; 2];                % local dofs of the "sides"
%         side_type = 'Point';            % sides are points
%         quad = ...                      % instance of Gauss quadrature class
%             mFEM.Gauss('Order', 1, 'Type', 'line');    
    end
    
    methods     
        function obj = Line2(id, nodes)
            %LINE2 Class constructor; calls base class constructor
            %
            % Syntax
            %   Line2(id, nodes)
            %
            % Description
            %   see mFEM.Element
            %
            % See Also mFEM.Element
            obj = obj@mFEM.cells.base.Cell(id, nodes); 
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
    
    methods (Static, Access = ?mFEM.Mesh)
        function [nodes,elements] = grid(x0,x1,xn)
            
        %     import mFEM.elements.base.* 
             import mFEM.elements.*
            
            %   type = mfilename('class');


            spmd  
                x = x0 : (x1-x0)/xn : x1;
                nodes = cell(length(x),1);
                elements = cell(length(x)-1,1);

                for i = 1:length(x);
                    nodes{i} = mFEM.elements.base.Node(i,x(i));

                    if i > 1;
                        elements{i-1} = Line2(i-1, nodes(i-1:i));
                    end
                end
                
                nodes = codistributed(nodes);
                elements = codistributed(elements);
             end
        end
    end
    
end