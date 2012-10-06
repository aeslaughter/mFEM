classdef Mesh < handle

    properties (SetAccess = private, GetAccess = public)
        n_elements = [];
        element = Quad4.empty;  % empty array of elements
        element_type = '';
        n_dim = [];        
        type = 'CG';
        space = 'scalar';
        map = struct('node', [], 'elem', uint32([]), 'dof', uint32([]));
        initialized = false;
    end
    
    properties (Access = private)
        CG_dof_map = uint32([]); % needed to compute neighbors for both CG and DG
    end
    
    methods (Access = public)
        function obj = Mesh(elem_name, varargin)
            
            % Assign element name
            obj.element_type = elem_name;

            % User provided the FE space ('scalar'(default) or 'vector')
            if nargin == 2;
                obj.space = varargin{1};
            end
            
            % User provided the FE type ('CG'(default) of 'DG')
            if nargin == 3;
                obj.type = varargin{2};
            end
            
        end
        
        function add_element(obj, varargin)
            obj.elements(end+1) = feval(obj.element_type, varargin{:});
        end
        
        function initialize(obj)
            obj.n_elements = length(obj.element); 
            obj.n_dim = obj.element(1).n_dim;
            obj.compute_dof_map();
            obj.find_neighbors();
            obj.initialized = true;
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
            % Calculates the global degree-of-freedom map
            % (this is used by both CG and DG for finding neighbors)
            [~,~,obj.CG_dof_map] = unique(obj.map.node, 'rows','stable');
            
            % Place the correct type of map in public property
            switch obj.type;
                case 'CG'; 
                    obj.map.dof = obj.CG_dof_map;
                case 'DG'; 
                     obj.map.dof = (1:size(obj.map.node,1))';
            end
            
            % Update the elements with the global dof
            for e = 1:obj.n_elements;
                elem = obj.element(e);
                elem.global_dof = obj.CG_dof_map(obj.map.elem == e,:);
            end
        end
        
        % This needs to be improved for working with 2D and 3D
        function plot(obj)
           
            figure; hold on;
            for e = 1:obj.n_elements;
                elem = obj.element(e);
                patch(elem.nodes(:,1), elem.nodes(:,2), 'b', 'FaceColor','none');
                text(mean(elem.nodes(1:2,1)),mean(elem.nodes(2:3,2)),num2str(e),...
                    'FontSize',14,...
                    'BackgroundColor','k','FontWeight','Bold',...
                    'HorizontalAlignment','center','Color','w');
            end
            
            N = length(unique(obj.map.dof));
            for i = 1:N;
                idx = obj.map.dof == i;
                x = unique(obj.map.node(idx,:),'rows');
               % plot(x(1),x(2),'ko');
                text(x(1),x(2),num2str(i),'FontSize',12,...
                    'BackgroundColor','b',...
                    'HorizontalAlignment','center','Color','w');
                 
            end
        end 
    end
    
    methods (Access = private)
        function gen2D_quad(obj, x0, x1, y0, y1, xn, yn, varargin) 

            % Generate the generic grid points
            x = x0 : (x1-x0)/xn : x1;
            y = y0 : (y1-y0)/yn : y1;

            % Loop through the grid, creating elements for each cell
            k = 0;
            %dof = 0;
            for j = 1:length(y)-1;           
                for i = 1:length(x)-1;
                    nodes(1,:) = [x(i), y(j)];
                    nodes(2,:) = [x(i+1), y(j)];
                    nodes(3,:) = [x(i+1), y(j+1)];
                    nodes(4,:) = [x(i), y(j+1)];
                   
                    k = k + 1;
                    obj.element(k) = feval(obj.element_type, k, nodes, obj.space);
                    nn = obj.element(k).n_nodes;
                    obj.map.elem(end+1:end+nn,:) = k;
                    obj.map.node(end+1:end+nn,:) = nodes;
                end
            end
            
            % Initialize
            obj.initialize();
        end 
                     
        function find_neighbors(obj)
            % Locates elements that share a side
            % (this could probably be written better)
            for e = 1:obj.n_elements;
                elem = obj.element(e);
                elem_dof  = elem.global_dof;
                
                for s = 1:elem.n_sides
                    % Determine side dof
                    side_dof = elem_dof(elem.side_nodes(s,:));
                    
                    % Collects index of neighbor elements (including corners)
                    idx = zeros(size(obj.CG_dof_map));
                    for i = 1:length(side_dof);
                        idx = idx + (obj.CG_dof_map == side_dof(i));
                    end
                    idx(obj.map.elem == e) = 0; %exclude current element

                    % Gathers neighbor element ids, only include as a neighbor
                    % if it shares multiple nodes (i.e., corners dont count)
                    x = obj.map.elem(logical(idx));
                    count = zeros(size(x));
                    for i = 1:length(x);
                        count(i) = sum(x(i) == x);
                    end
                    id = unique(x(count > 1));
                    
                    % Insert neigbhor elements into current element
                    neighbors = feval([obj.element_type,'.empty']);
                    for i = 1:length(id);
                       neighbors(i) = obj.element(id(i));
                    end
                    elem.neighbor(s).element = neighbors;   
                end
            end
        end  
    end
end