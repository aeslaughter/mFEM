classdef Quad4 < mFEM.cells.base.Cell
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
    properties (Access = protected) 
        n_sides = 4;                        % no. of sides
        side_ids = [1,2; 2,3; 3,4; 4,1];    % define the side dofs 
%         side_type = 'Line2';                % side is 2-node line element
%         quad = ...                          % quadrature rules
%             mFEM.Gauss('order', 2, 'type', 'quad');       
     end
    
    % Define the Quad4 constructor
    methods
        function obj = Quad4(varargin)
           % Class constructor; calls base class constructor
           obj = obj@mFEM.cells.base.Cell(varargin{:}); 
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
    
    methods (Static, Access = ?mFEM.Mesh)
        function [nodes,elements] = grid(x0,x1,y0,y1,xn,yn)

            spmd
                id = 0;
                n_nodes = 0;
                N = zeros(4,xn*yn);

                for i = 1:xn;
                    for j = 1:yn;
                        id = id + 1;
                        if i == 1 && j == 1;
                            N(:,id) = n_nodes+1:n_nodes+4;
                            n_nodes = n_nodes + 4;
                        elseif i == 1 && j > 1;
                            b = id-1; 
                            N(:,id) = [N(b,4), N(3,b), n_nodes+1, n_nodes+2];
                            n_nodes = n_nodes + 2;
                        elseif i > 1 && j == 1;
                            l = id-yn;
                            N(:,id) = [N(2,l), n_nodes+1, n_nodes+2, N(3,l)];
                            n_nodes = n_nodes + 2;
                        else
                            b = id-1;
                            l = id-yn;
                            N(:,id) = [N(2,l), N(3,b), n_nodes+1, N(3,l)];
                            n_nodes = n_nodes + 1;
                        end
                    end
                end
                
                N = codistributed(N);
            end
              
            spmd
                X = codistributed(x0 : (x1-x0)/xn : x1);
                Y = codistributed(y0 : (y1-y0)/yn : y1);
               
                x = getLocalPart(X);
                y = getLocalPart(Y);
                

                local = cell(length(x)*length(y),1);
                k = 0;
                for i = 1:length(x);
                    for j = 1:length(y);
                        k = k + 1;
%                         local{k} = mFEM.elements.base.Node([x(i),y(j)]);
                    end
                end
                
                codist = codistributor1d(1,codistributor1d.unsetPartition, [n_nodes,1]);
                nodes = codistributed.build(local,codist);
            end

            
            
%           nodes = {};
          elements = {};
            
%           spmd
%               N_local = getLocalPart(N);
%               id_local = globalIndices(N,2);
%               codist = getCodistributor(N);
%               id = cumsum(codist.Partition);
% 
%                 for i = 1:xn;
%                     id = id + 1;
% 
%                     if i == 1 && j == 1;
%                        nodes{n+1} = mFEM.elements.base.Node(n+1,[x(i),y(j)]);
%                        nodes{n+2} = mFEM.elements.base.Node(n+2,[x(i+1), y(j)]);
%                        nodes{n+3} = mFEM.elements.base.Node(n+3,[x(i+1), y(j+1)]);
%                        nodes{n+4} = mFEM.elements.base.Node(n+4,[x(i), y(j+1)]);             
%                         N(id,:) = n+1:n+4;
%                         n = n + 4;
%                     elseif i == 1 && j > 1;
%                         b = id-1;
%                        nodes(n+1) = {mFEM.elements.base.Node(n+1,[x(i+1), y(j+1)])};
%                        nodes(n+2) = {mFEM.elements.base.Node(n+2,[x(i), y(j+1)])};  
%                         N(id,:) = [N(b,4), N(b,3), n+1, n+2];
%                         n = n + 2;
%                     elseif i > 1 && j == 1;
%                         l = id-yn;
%                        nodes(n+1) = {mFEM.elements.base.Node(n+1,[x(i+1), y(j)])};
%                        nodes(n+2) = {mFEM.elements.base.Node(n+2,[x(i+1), y(j+1)])};
%                         N(id,:) = [N(l,2), n+1, n+2, N(l,3)];
%                         n = n + 2;
%                     else
%                         b = id-1;
%                         l = id-yn;
%                         nodes(n+1) = {mFEM.elements.base.Node(n+1,[x(i+1), y(j+1)])}; 
%                         N(id,:) = [N(l,2), N(b,3), n+1, N(l,3)];
%                         n = n + 1;
%                     end
% 
%                     nodes = codistributed(nodes);
%                     elements{id} = mFEM.elements.Quad4(id,nodes(N(id,:)));
%                 end
%           end
          
          
          
          
          
          
%             x = x0 : (x1-x0)/xn : x1;
%             y = y0 : (y1-y0)/yn : y1;
% 
% %             spmd
%             
%                 % Loop through the grid, creating elements for each cell
%                 id = 0;
%                 n = 0;
% 
%                 N = zeros(xn*yn,4);
%                 nodes = cell(length(x)*length(y),1); 
%                 elements = cell(xn*yn,1);
% 
%                 for i = 1:xn;
%                     for j = 1:yn;
%                         id = id + 1;
%                         
%                         if i == 1 && j == 1;
%                            % nodes{n+1} = mFEM.elements.base.Node(n+1,[x(i),y(j)]);
%                            % nodes{n+2} = mFEM.elements.base.Node(n+2,[x(i+1), y(j)]);
%                            % nodes{n+3} = mFEM.elements.base.Node(n+3,[x(i+1), y(j+1)]);
%                            % nodes{n+4} = mFEM.elements.base.Node(n+4,[x(i), y(j+1)]);             
%                             N(id,:) = n+1:n+4;
%                             n = n + 4;
%                         elseif i == 1 && j > 1;
%                             b = id-1;
%                            % nodes(n+1) = {mFEM.elements.base.Node(n+1,[x(i+1), y(j+1)])};
%                            % nodes(n+2) = {mFEM.elements.base.Node(n+2,[x(i), y(j+1)])};  
%                             N(id,:) = [N(b,4), N(b,3), n+1, n+2];
%                             n = n + 2;
%                         elseif i > 1 && j == 1;
%                             l = id-yn;
%                            % nodes(n+1) = {mFEM.elements.base.Node(n+1,[x(i+1), y(j)])};
%                            % nodes(n+2) = {mFEM.elements.base.Node(n+2,[x(i+1), y(j+1)])};
%                             N(id,:) = [N(l,2), n+1, n+2, N(l,3)];
%                             n = n + 2;
%                         else
%                             b = id-1;
%                             l = id-yn;
%                             %nodes(n+1) = {mFEM.elements.base.Node(n+1,[x(i+1), y(j+1)])}; 
%                             N(id,:) = [N(l,2), N(b,3), n+1, N(l,3)];
%                             n = n + 1;
%                         end
% 
%                         %nodes = codistributed(nodes);
%                         %elements{id} = mFEM.elements.Quad4(id,nodes(N(id,:)));
%                     end
%                 end
%             end
        end
    end    
end