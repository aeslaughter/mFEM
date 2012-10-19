classdef Tri6 < mFEM.Element
    %Tri6 6-node triangle element
    % 
    %        xi2  
    %         ^
    %         |
    %               
    %         3
    %         | \  
    %         |  \
    %         |   \
    %     (3) 6    5  (2)      
    %         |     \           
    %         |      \    
    %         1---4---2  ---> xi1
    %            (1) 
    %
    % http://mmc.geofisica.unam.mx/Bibliografia/Matematicas/EDP/MetodosNumericos/FEM/IFEM.Ch24.pdf

    % Define the inherited abstract properties
    properties (SetAccess = protected, GetAccess = public)
        n_sides = 3;                % no. of sides
        lims = [0,1];               % limits of xi1 and xi2        
        side_dof = [1,4,2; 2,5,3; 3,6,1]; % define the side dofs 
        side_type = 'Linear2';      % uses 2-node linear element for sides
    end
    
    % Define the Quad4 constructor
    methods
        function obj = Tri6(id, nodes, varargin)
           % Class constructor; calls base class constructor
           
           % Test that nodes is sized correctly
           if ~all(size(nodes) == [3,2]) && ~all(size(nodes) == [6,2]);
                error('Tri6:Tri6','Nodes not specified correctly; expected a [3x2] or [6x2] array, but recieved a [%dx%d] array.', size(nodes,1), size(nodes,2));
           end
           
           % Special case when only nodes 1,2,3 are given
           if all(size(nodes) == [3,2])
               nodes(4,:) = mean(nodes(1:2,:));
               nodes(5,:) = mean(nodes(2:3,:));
               nodes(6,:) = mesh(nodes([3,1],:));
           end

           % Call the base class constructor
           obj = obj@mFEM.Element(id, nodes, varargin{:});
        end
    end
    
    % Define the inherited abstract methods (protected)
    methods (Access = protected)    
        
        function J = jacobian(obj, xi1, xi2)
            % Define the Jacobain
            x = obj.nodes(:,1);
            y = obj.nodes(:,2);
            xi3 = 1 - xi1 - xi2;
            
            Jx1 = x(1)*(4*xi1-1) + 4*(x(4)*xi2 + x(6)*xi3);
            Jx2 = x(2)*(4*xi2-1) + 4*(x(5)*xi3 + x(4)*xi1);
            Jx3 = x(3)*(4*xi3-1) + 4*(x(6)*xi1 + x(5)*xi2);
            Jy1 = y(1)*(4*xi1-1) + 4*(y(4)*xi2 + y(6)*xi3);
            Jy2 = y(2)*(4*xi2-1) + 4*(y(5)*xi3 + y(4)*xi1);
            Jy3 = y(3)*(4*xi3-1) + 4*(y(6)*xi1 + y(5)*xi2);            
            
            J = [1,1,1; Jx1, Jx2, Jx3; Jy1, Jy2, Jy3];
        end
        
        function N = basis(~, xi1, xi2)
            % Returns a row vector of local shape functions     
            xi3 = 1 - xi1 - xi2;
            
            N(1) = xi3*(2*xi3-1);
            N(2) = xi1*(2*xi1-1);
            N(3) = xi2*(2*xi2-1);
            N(4) = 4*xi3*xi1;
            N(5) = 4*xi1*xi2;
            N(6) = 4*xi3*xi2;
          end

        function B = grad_basis(obj, xi1, xi2) 
            % Gradient of shape functions    
            xi3 = 1 - xi1 - xi2;
            x = obj.nodes(:,1);
            y = obj.nodes(:,2);
            
            Dx4 = x(4) - 1/2*(x(1) + x(2));
            Dy4 = y(4) - 1/2*(y(1) + y(2));
            Dx5 = x(5) - 1/2*(x(2) + x(3));
            Dy5 = y(5) - 1/2*(y(2) + y(3));           
            Dx6 = x(6) - 1/2*(x(3) + x(1));
            Dy6 = y(6) - 1/2*(y(3) + y(1));                
            
            xx = @(i,j) obj.nodes(i,1) - obj.nodes(j,1);
            yy = @(i,j) obj.nodes(i,2) - obj.nodes(j,2);
            
            Jx21 = xx(2,1) + 4*(Dx4*(xi1-xi2) + (Dx5-Dx6)*xi3);
            Jx32 = xx(3,2) + 4*(Dx5*(xi2-xi3) + (Dx6-Dx4)*xi1);
            Jx13 = xx(1,3) + 4*(Dx6*(xi3-xi1) + (Dx4-Dx5)*xi2);
            Jy12 = yy(1,2) + 4*(Dy4*(xi2-xi1) + (Dy6-Dy5)*xi3);
            Jy23 = yy(2,3) + 4*(Dy5*(xi3-xi2) + (Dy4-Dy6)*xi1);
            Jy31 = yy(3,1) + 4*(Dy6*(xi1-xi3) + (Dy5-Dy4)*xi2);
            
            B(:,1) = [(4*xi1-1)*Jy23, (4*xi1-1)*Jx32];
            B(:,2) = [(4*xi2-1)*Jy31, (4*xi2-1)*Jx13];   
            B(:,3) = [(4*xi3-1)*Jy12, (4*xi3-1)*Jx21];
            B(:,4) = [4*(xi2*Jy23+xi1*Jy31), 4*(xi2*Jx32+xi1*Jx13)];
            B(:,5) = [4*(xi3*Jy31+xi2*Jy12), 4*(xi3*Jx13+xi2*Jx21)]; 
            B(:,6) = [4*(xi1*Jy12+xi3*Jy23), 4*(xi1*Jx21+xi3*Jx32)];              
        end
        
        function GN = local_grad_basis(~, varargin)
            error('Tri6:local_grad_basis', 'Function not defined for the %s element, the B matrix is computed directly.', class(obj));
        end
    end
end