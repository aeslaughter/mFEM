classdef Mesh < handle

    properties (SetAccess = private, GetAccess = public)
        n_elements = [];
        element_type = '';
        n_dim = [];        
        type = 'CG';
        node_map = [];
        n_nodes = [];
        dof_map = uint32([]);
        initialized = false;
    end
    

    properties (Access = private)
        elements_ = {};
    end
    
    methods (Access = public)
        function obj = Mesh(elem_name, varargin)
            
            if nargin >= 1;
                obj.element_type = elem_name;
                obj.n_dim = 2;
            end
            
            if nargin == 2
                obj.type = varargin{1};
            end
            
        end
        
        function initialize(obj)
           obj.compute_dof_map();
           % obj.find_neighbors();
           obj.initialized = true;
        end
        
        function e = element(obj, id)
            e = obj.elements_{id};
        end
        
        
        function gen2D(obj, x0, x1, y0, y1, xn, yn, varargin)  
            
            switch obj.element_type;
                case 'Quad4';
                    obj.gen2D_quad(x0, x1, y0, y1, xn, yn, varargin{:});
                otherwise
                    error('Mesh.gen2D is not supported for this element');
            end
        end
       
        function compute_dof_map(obj)
           
            switch obj.type;
                case 'CG'; obj.CG_compute_dof_map();
                case 'DG'; obj.DG_compute_dof_map();
            end
        end
        
        
    end
    
    methods (Access = private)
        function gen2D_quad(obj, x0, x1, y0, y1, xn, yn, varargin) 

            % Generate the generic grid
            xn = (x1 - x0) / xn;
            yn = (y1 - y0) / yn;
            [X,Y] = meshgrid(x0:xn:x1, y0:yn:y1);
            [m,n] = size(X);
        
            % Loop through the grid, creating elements for each cell
            k = 0;
            %dof = 0;
            for i = 1:m-1;
                for j = 1:n-1;
                    x = [X(i,j), X(i,j+1), X(i+1, j+1), X(i+1,j)];
                    y = [Y(i,j), Y(i,j+1), Y(i+1, j+1), Y(i+1,j)];
                    [x,y] = obj.cell2nodes2D(x,y);
                    n = length(x); % no. of nodes on element
                    
                    for e = 1:size(x,1);
                 
                        k = k + 1;
                        obj.elements_{k} = feval(obj.element_type, k, x(e,:), y(e,:));
                        
                        obj.node_map(end+1:end+n,:) = [x', y'];

%                         n_dof = elem.n_shape*elem.n_dof;
%                         elem.global_dof = dof+1:dof+n_dof;
%                         dof = dof + n_dof;      
%                         obj.elements_{k} = elem;                           
                    end
                end
            end
            
            % Define no. of elements
            obj.n_elements = k;
            
            % Compute dof map
            obj.compute_dof_map;
           % obj.initialized = true;  

        end 
        
        function [nx,ny] = cell2nodes2D(obj,x,y)
            % Converts x,y grid (4 points) to nodal values
            
            switch obj.element_type;
                case 'Quad4';
                    nx = x;
                    ny = y;
%                case 'Quad8';
%                     nx = [x, mean(x(1:2)), mean(x(2:3)), mean(x(3:4)), mean(x([1,4]))];
%                     ny = [y, mean(y(1:2)), mean(y(2:3)), mean(y(3:4)), mean(y([1,4]))];
            end
            
        end   
         
        function CG_compute_dof_map(obj)
            % Computes the global dof map for continous elements
            
            % Identify the unique nodes
            [C,~,~] = unique(obj.node_map,'rows');
            n = size(C,1); % no. of unique nodes
            
            % Determine the number of rows and cols in map
            [mr,mc] = size(obj.node_map);

            % Loop through unique nodes
            for i = 1:n;
                
                % Matchs x, y, z values for each unique node
                idx = zeros(mr,1);
                for j = 1:mc;
                    idx = idx + (obj.node_map(:,j) == C(i,j));
                end

                % Assigns node no. to dof_map vector
                obj.dof_map(idx == mc,1) = uint32(i);
            end
            
            % Assigns the total no. of nodes to the object property
            obj.n_nodes = n;
        end
        
        function DG_compute_dof_map(obj)
            % Computes global dof map for discontinous elements

            % Each node is treated independantly
            obj.dof_map = (1:size(obj.node_map,1))';
        end
        
  
    end
end