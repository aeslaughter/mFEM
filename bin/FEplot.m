%% FEplot
% Creates graphs for finite element data derived from FEMESH.
%
%% Syntax
%  FEplot(obj)
%  FEplot(obj, data)
%  FEplot(...,'PropertyName',<PropertyValue>)
%
%% Descrtiption
%
%% FEplot Property Descriptions
%
%  ShowNodes
%       {true} | false
%
%       Toggles the display of plotted dots at the nodes
%
%  ElementLabels
%       true (default w/o data) | false (default w/ data)
%
%       Toggles the appearence of element number labels.
%       
%   NodeLabels
%       true (default w/o data) | false (default w/ data)
%
%       Toggles the appearence of node number labels.
%
%   NewFigure
%       true (default w/o data) | false (default w/ data)
%
%       Toggles the creation of a new figure, if it is set to false then
%       gcf is used, unless the handle is specified via FigureHandle
%
%   AxesHandle
%       axes handle
%
%       Add the plot to the user specified axes handle. 
%
%   FigureHandle
%       figure handle
%
%       Add the plot to the user specified figure handle. Note, NewFigure
%       takes presidence over FigureHandle.
%
%% See also
% FEmesh

function FEplot(obj, varargin)

    % Check that mesh is initialized
    if (~obj.initialized)
        error('ERROR: FEmesh must be initialized before graphs may be plotted.');
    end
   
    % Collect the input 
    opt = parse_input(obj, varargin{:});

    % Plot the data
    build_plot(obj, opt);
    
    % Add the element labels
    add_element_labels(obj, opt);
       
    % Add the nodel lables
    add_node_labels(obj, opt)
    
    % Resize the figure
    box on;
end 

function opt = parse_input(obj, varargin)
    %PARSE_INPUT Function for parsing input data to plot function
    %
    % See FEPLOT

    % Define user properties
    opt.data = []; 
    opt.shownodes = false;
    opt.elementlabels = true;
    opt.nodelabels = true;
    opt.newfigure = true;
    opt.figurehandle = handle.empty;
    opt.axeshandle = handle.empty; 
    
    % The first input is always the data
    if nargin > 1 && isnumeric(varargin{1});
        opt.data =  varargin{1};
        
        % Check that the data is sized correctly
        if ~isempty(opt.data) && length(opt.data) ~= obj.n_dof;
            error('FEmesh:plot', 'Data not formated correctly, it must be a vector of numbers of length %d.', obj.n_dof);
        end
    end

    % Parse the input
    if nargin >= 2;

        % When using data, disable the labels by default and do not
        % create a new figure with each call
        opt.elementlabels = false;
        opt.nodelabels = false;
        opt.newfigure = false;

        % Collect options supplied by the user
        opt = gather_user_options(opt, varargin{2:end});
    end
    
    % Set/Create the figure handle
    if opt.newfigure;
        figure('color','w'); hold on;   
    elseif ishandle(opt.figurehandle)
        figure(opt.figurehandle,'color','w'); hold on;
    end
    
    % Set the axes handle
    if ishandle(opt.axeshandle);
        axes(opt.axeshandle);
    end
end

function build_plot(obj,opt)
    %BUILD_PLOT

    % 1D case (data is added via cell_node_data)
    if obj.n_dim == 1;
        cnode = cell_node_data(obj, opt);
        h = plot(cnode{:}); 
    
    % 2D case    
    else
        h = zeros(obj.n_elements,1);
        for e = 1:obj.n_elements;
            % Current element
            elem = obj.element(e);
            
            % Get the nodal data
            cnode = cell_node_data(elem, opt);
        
            if ~isempty(opt.data); % with data
                h(e) = patch(cnode{:}, opt.data(elem.get_dof));
            else % mesh only
                h(e) = patch(cnode{:}, 'k', 'FaceColor','none');
            end
        end
    end

    % Show the actual nodes
    if opt.shownodes
        set(h, 'Marker','o', 'MarkerSize', 8, 'MarkerEdgeColor', 'auto',...
            'MarkerFaceColor', 'auto');
    end
end
    
function varargout = cell_node_data(obj, opt)
    %GET_NODE_DATA Creates column cell arrays of nodal data

    % Get the node data (all elements; FEmesh object)
    if isa(obj,'mFEM.FEmesh');
        node = obj.map.node;
        
        % Adjust node data for CG or DG elements
        if strcmpi(obj.opt.type,'CG');
            node = unique(node,'rows','R2012a');
        end
        
    % Get the node data (single element; Element object)    
    elseif isa(obj,'mFEM.Element');
        node = obj.nodes;
        
        if ~isempty(obj.node_plot_order);
           node = node(obj.node_plot_order,:);
        end
        
    end
    
    % Create y-direction for 1D case
    if obj.n_dim == 1;  
        % If data is given, use this for the y position
        if isempty(opt.data); 
            node(:,2) = zeros(size(node(:,1)));
        else
            node(:,2) = opt.data(obj.get_dof);
        end
    end

    % Create cell structure of nodes
    varargout{1} = num2cell(node,1); % col. cell of nodes
    
    % Compute  mean locations of the nodes (only for element case)
    if isa(obj,'mFEM.Element') && nargout == 2;
    	varargout{2} = num2cell(mean(node,1)); % location of nodes
    end
end

function add_element_labels(obj, opt)
    %ADD_ELEMENT_LABEL
    
    % Only continue if the labels are wanted
    if ~opt.elementlabels; return; end
    
    % Loop through all of the elements
    for e = 1:obj.n_elements;

        % The current element
        elem = obj.element(e);

        % Get the label position
        [~,cntr] = cell_node_data(elem, opt);

        % Add the label
        text(cntr{:}, num2str(e),'FontSize',14,...
            'BackgroundColor','k','FontWeight','Bold',...
            'Color','w','HorizontalAlignment','center');
    end
end

function add_node_labels(obj, opt)
    %ADD_NODE_LABELS
    
    % Return if the node labels are not desired
    if ~opt.nodelabels;
        return;
    end

    % The node labels are placed differently for CG and DG meshes
    switch obj.opt.type;
        
        % Continous finite element labels
        case 'CG'; 

            % Extract the unique nodes
            N = unique(obj.map.node,'rows','stable');
            
            % Adjust for the 1D case
            if obj.n_dim == 1;
                
                % If data is given, use this for the y position
                if isempty(opt.data); 
                    N(:,2) = zeros(size(N(:,1)));
                else
                    N(:,2) = opt.data;
                end
            end
            
            % Loop through the unique nodes and create the labels=
            for i = 1:length(N);
                node = num2cell(N(i,:));
                text(node{:},num2str(i),'FontSize',12,'Color','w',...
                    'BackgroundColor','b','HorizontalAlignment','center');
            end
            
        % Discontinous finite element labels
        case 'DG'
            
            % Loop through the elements
            for e = 1:obj.n_elements;
                
                % Get the current elem
                elem = obj.element(e);
                
                % Get the nodes locations
                node = cell_node_data(elem, opt);
                N = cell2mat(node);
                  
                % Get the global degrees of freedom for the element
                dof = elem.get_dof();
                
                % Loop through the nodes and add the text
                for i = 1:size(N,1);
                    
                    % Locate the nodal labels
                    [X,Y,n] = locate_node_labels(elem, N(i,:));
                    
                    % Create label
                    node = num2cell(n);
                    text(node{:},num2str(dof(i)),'FontSize',12,'Color','w',...
                        'BackgroundColor','b','HorizontalAlignment',X,...
                        'VerticalAlignment',Y);
                end  
            end     
    end
end

function [X,Y,N] = locate_node_labels(elem, N)
    %LOCATE_NODE_LABELS
       
    % Get the limiting positions 
    if elem.n_dim >= 1;
        
        
        % Get the real locations of the far extents of element
        [x0] = elem.get_position(elem.lims(1),elem.lims(1));
        [x1] = elem.get_position(elem.lims(2),elem.lims(1));
        padX = 0.05*(x1-x0);

        % Set location of label (x-direction)
        if N(1) == x0; 
            X = 'left';
            N(1) = N(1) + padX;
        elseif N(1) == x1; 
            X = 'right';
            N(1) = N(1) - padX;
        else X = 'center';
        end
        
        % Set the location of Y
        Y = 'middle';
    end
    
    if elem.n_dim >= 2;
        % Get the real locations of the far extents of element
        [~,y0] = elem.get_position(elem.lims(1),elem.lims(1));
        [~,y1] = elem.get_position(elem.lims(2),elem.lims(2));
        padY = 0.05*(y1-y0);    
        
        % Set location of label (x-direction)
        if N(2) == y0; 
            Y = 'top';
            N(2) = N(2) - padY;
        elseif N(2) == x1; 
            Y = 'bottom';
            N(2) = N(2) + padX;
        else Y = 'middle';
        end
    end

end
