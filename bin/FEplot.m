% FEPLOT Creates graphs for finite element data derived from FEMESH.
%
% Syntax
%  FEplot(obj)
%  FEplot(obj, data)
%  FEplot(...,'PropertyName',<PropertyValue>)
%
% Descrtiption
%   FEplot(obj) creates a plot of the mesh object (obj) and includes that
%   includes the element and degree of freedom labels.
%
%   FEplot(obj, data) creates a plot for the mesh object (obj) and the data 
%   supplied, the data must be a column vector of with the number of values 
%   equal to the number of degrees of freedom for the mesh.
%
%   FEplot(...,'PropertyName',PropertyValue,...) allows for the
%   plot created to be customised, for a list of properties below. Note 
%   that all true/false may be input using a flag instead of the pairing. 
%   For example the following two calls to FEplot are the same:
%       FEplot(obj,data,'-Deform');
%       FEplot(obj,data,'Deform',true);
%   This flag style input simply reverses the state from the default.
%
% FEplot Property Descriptions
%  Deform
%       true | {false}
%       Toggles the use to treat the input data as deformation
%
%  Component
%       interger
%       Gives the vector component to plot as the countour variable, for
%       example 1 would used the x-component of a vector valued supplied
%       data. If the value is empty (the default) the magnitude is shown,
%       if it is set to NaN zeros are used. The deformation property shows 
%       the deformed structure visually.
%
%  ShowNodes
%       {true} | false
%       Toggles the display of plotted dots at the nodes
%
%  ElementLabels
%       true (default w/o data) | false (default w/ data)
%       Toggles the appearence of element number labels.
%       
%   NodeLabels
%       true (default w/o data) | false (default w/ data)
%       Toggles the appearence of node number labels.
%
%   New
%       true (default w/o data) | false (default w/ data)
%       Toggles the creation of a new figure, if it is set to false then
%       gcf is used, unless the handle is specified via FigureHandle
%
%   Axes
%       axes handle
%       Add the plot to the user specified axes handle. 
%
%   Figure
%       figure handle
%       Add the plot to the user specified figure handle. Note, NewFigure
%       takes presidence over FigureHandle.
%
%   Plot
%       cell array
%       Use this to pass commands directly to the plot (1D), for example
%           FEplot(...,'Plot',{'Xlim',[0,1]});
%       executes the following:
%           plot(...,'Xlim',[0,1]);
%       This is the last command to be applied to the plot, thus it will
%       override other plot related settings.
%
%   Patch
%       cell array
%       Use this to pass commands directly to patch (2D & 3D), for example
%           FEplot(...,'Patch',{'FaceColor','k'});
%       executes the following:
%           patch(...,'FaceColor','k');
%       This is the last command to be applied to the patch, thus it will
%       override other patch related settings.
%
% See also
% FEmesh
%
%----------------------------------------------------------------------
% Copyright 2012 Andrew E. Slaughter
% This software is for educational purposes only and may not be used
% without written permession.
%----------------------------------------------------------------------

function FEplot(obj, varargin)

    % Check that mesh is initialized
    if (~obj.initialized)
        error('ERROR: FEmesh must be initialized before graphs may be plotted.');
    end
   
    % Collect the input 
    opt = parse_input(obj, varargin{:});

    % Plot the data according the spacial dimensions
    if obj.n_dim == 1 && obj.n_dof_node == 1;
        h = plot1D_scalar(obj, opt);
        
    elseif obj.n_dim == 2 && obj.n_dof_node == 1;
        h = plot2D_scalar(obj, opt);
        
    elseif obj.n_dim == 2 && obj.n_dof_node == 2;
        h = plot2D_vector(obj, opt);
    end
    
    apply_plot_options(h, obj, opt);
    %build_plot(obj, opt);
    
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
    opt.deform = false;
    opt.scale = 'auto';
    opt.component = [];
    opt.polar = false;
    opt.colorbar = '';
    opt.shownodes = false;
    opt.elementlabels = true;
    opt.nodelabels = true;
    opt.new = true;
    opt.figure = handle.empty;
    opt.axes = handle.empty; 
    opt.plot = {};
    opt.patch = {};
    opt.hold = true;
    
    % Account for first input containing the data
    if nargin > 1 && isnumeric(varargin{1});
        opt.data =  varargin{1};
        start_idx = 2;
        
        % Check that the data is sized correctly
        if ~isempty(opt.data) && length(opt.data) ~= obj.n_dof;
            error('FEmesh:plot', 'Data not formated correctly, it must be a vector of numbers of length %d.', obj.n_dof);
        end
    else
        start_idx = 1;
    end

    % When using data, disable the labels by default and do not
    % create a new figure with each call
    if ~isempty(opt.data);
        opt.elementlabels = false;
        opt.nodelabels = false;
        opt.newfigure = false;
    end
    
    % Collect options supplied by the user
    opt = gather_user_options(opt, varargin{start_idx:end}); 

    % Set/Create the figure handle
    if opt.newfigure;
        figure('color','w'); hold on;   
    elseif ishandle(opt.figure)
        figure(opt.figure,'color','w'); hold on;
    end
    
    % Set the axes handle
    if ishandle(opt.axes);
        axes(opt.axes);
    end
    
    % Determine the scale factor
    if ~isempty(opt.data) && opt.deform && strcmpi(opt.scale,'auto');
        % Determine the scaling
        n = max(max(floor(log10(abs(opt.data)))));
        opt.scale = 10^(-sign(n)*(abs(n)-1));   
    end
end

function h = plot1D_scalar(obj, opt)
    %PLOT1D create a 1D plot

    % Initialize the handle output
    X = [];
    Y = [];
    
    % Generate the graph
    for e = 1:obj.n_elements;

        % The current element
        elem = obj.element(e);
        
        % The node positions
        x = elem.nodes;
   
        % Gather y-axis data
        if ~isempty(opt.data);        
            dof = elem.get_dof();
            y = opt.data(dof);
        else
            y = zeros(size(x));
        end
        
        % Get the plotting order
        order = elem.node_plot_order;
        if ~isempty(order);
            order = elem.node_plot_order;
            x = x(order);
            y = y(order);
        end

        % Append the output
        X = [X;x];
        Y = [Y;y];
    end

    % Create the graph
    h = plot(X,Y); hold on;
end

function h = plot2D_scalar(obj, opt)
    %PLOT1D create a 2D plot

    % Initialize the x and y values
    h = zeros(obj.n_elements,1);
    
    % Generate the graph
    for e = 1:obj.n_elements;

        % The current element
        elem = obj.element(e);
        
        % The node positions
        x = elem.nodes(:,1);
        y = elem.nodes(:,2);
        
        % Gather y-axis data
        z = [];
        if ~isempty(opt.data)
            dof = elem.get_dof();
            z(:,1) = opt.data(dof);
        end
        
        % Apply the plotting order
        order = elem.node_plot_order;
        if ~isempty(order);
            x = x(order);
            y = y(order);
            if ~isempty(z);
                z = z(order);
            end
        end

        % Create the graph
        if isempty(z);
            h = patch(x, y, 'w');
        else
            h = patch(x, y, z);
        end
    end
end

function h = plot2D_vector(obj, opt)
    %PLOT1D create a 2D plot for vector spaces

    % Initialize the x and y values
    h = zeros(obj.n_elements,1);
    
    % Generate the graph
    for e = 1:obj.n_elements;

        % The current element
        elem = obj.element(e);
        
        % The node positions
        x = elem.nodes(:,1);
        y = elem.nodes(:,2);
        
        % Adjust for polar coordinates
        if opt.polar;
            [x,y] = pol2cart(x,y);
        end

        % Gather y-axis data
        z = [];
        if ~isempty(opt.data)
            dof = elem.get_dof();
            z = opt.data(dof);
            zz(:,1) = z(1:2:end);
            zz(:,2) = z(2:2:end);           
            
            if isnan(opt.component);
                z = nan(size(x));
            elseif isempty(opt.component);
                z = sqrt(zz(:,1).^2 + zz(:,2).^2);
            elseif isnumeric(opt.component) && isscalar(opt.component);
                z = zz(:,opt.component);
            else
                error('FEplot:plot2D_vector','Error with Component property value');
            end
            
             % Adjust nodal position if deformed shape is desired
            if opt.deform
                zx = zz(:,1); zy = zz(:,2);
                
                % Adjust for polar cordinates on z
                if opt.polar;
                    [zx, zy] = pol2cart(zx,zy);
                end
                
                % Adjust the nodal values
                x = x + opt.scale*zx;
                y = y + opt.scale*zy;
            end
        end

        % Apply the plotting order
        order = elem.node_plot_order;
        if ~isempty(order);
            x = x(order);
            y = y(order);
            if ~isempty(z);
                z = z(order);
            end
        end

        % Create the graph
        if isempty(z);
            h(e) = patch(x, y, 'w');
        else
            h(e) = patch(x, y, z);
        end
        
        if ~isempty(opt.colorbar) && ischar(opt.colorbar);
            cbar = colorbar;
            set(get(cbar,'YLabel'),'String',opt.colorbar);
        end
    end
end

function apply_plot_options(h, obj, opt)
    % APPLY_PLOT_OPTIONS
    
    % Show the nodes as empty circles
    if opt.shownodes
        set(h, 'Marker','o', 'MarkerSize', 8, 'MarkerEdgeColor', 'auto',...
            'MarkerFaceColor', 'auto');
    end
    
    % Hold the plot
    if opt.hold;
        hold on; 
    else
        hold off; 
    end
   
    % Add the custom plot options
    if ~isempty(opt.plot) && obj.n_dim == 1;
        set(h, opt.plot{:}); 
    end
    
    % Add the custom patch options
    if ~isempty(opt.patch) && obj.n_dim > 1;
        set(h, opt.patch{:}); 
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
        nodes = elem.nodes;

        % Adjust for polar coordinates
        if opt.polar;
            [nodes(:,1),nodes(:,2)] = pol2cart(nodes(:,1),nodes(:,2));
        end
        
        % Determine the center position
        cntr = num2cell(mean(nodes,1));

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
            nodes = unique(obj.map.node,'rows','stable');
            
            % Adjust for polar coordinates
            if opt.polar;
                [nodes(:,1),nodes(:,2)] = pol2cart(nodes(:,1),nodes(:,2));
            end
            
            % Adjust for the 1D case
            if obj.n_dim == 1;
                
                % If data is given, use this for the y position
                if isempty(opt.data); 
                    nodes(:,2) = zeros(size(nodes(:,1)));
                else
                    nodes(:,2) = opt.data;
                end
            end
            
            % Loop through the unique nodes and create the labels=
            for i = 1:length(nodes);
                node = num2cell(nodes(i,:));
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
                nodes = cell2mat(node);
                  
                % Get the global degrees of freedom for the element
                dof = elem.get_dof();
                
                % Loop through the nodes and add the text
                for i = 1:size(nodes,1);
                    
                    % Locate the nodal labels
                    [X,Y,n] = locate_node_labels(elem, nodes(i,:));
                    
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
