% FEPLOT Creates graphs for finite element data derived from FEMESH.
%
% Syntax
%  FEplot(obj, [])
%  FEplot(obj, data)
%  FEplot(...,'PropertyName',<PropertyValue>)
%
% Descrtiption
%   FEplot(obj, []) creates a plot of the mesh object (obj) and includes that
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
%  Scale
%       {'auto'} | numeric
%       Set the deformation scale factor, which is only appliciable if the
%       deformation property is used. My default a factor is automatically
%       computed based on the order of the data. It may be changed by
%       specifing a value (e.g., FEplot(obj,'Scale',100)). To see the
%       actual displacement set the value to 1.
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
%  NodeLabels
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

function FEplot(obj, data, varargin)

    % Check that mesh is initialized
    if (~obj.initialized)
        error('ERROR: FEmesh must be initialized before graphs may be plotted.');
    end
   
    % Collect the input 
    opt = parse_input(obj, data, varargin{:});

    % Plot the data according the spacial dimensions
    if obj.n_dim == 1 && obj.n_dof_node == 1;
        h = plot1D_scalar(obj, opt);
        
    elseif obj.n_dim == 1 && obj.n_dof_node == 2;
        h = plot1D_vector(obj, opt);
        
    elseif obj.n_dim == 2 && obj.n_dof_node == 1;
        h = plot2D_scalar(obj, opt);
        
    elseif obj.n_dim == 2 && obj.n_dof_node == 2;
        h = plot2D_vector(obj, opt);
        
    elseif obj.n_dim == 3 && obj.n_dof_node == 1;
        h = plot3D_scalar(obj, opt);
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

function opt = parse_input(obj, data, varargin)
    %PARSE_INPUT Function for parsing input data to plot function
    %
    % See FEPLOT

    % Define user properties
    opt.data = data; 
    opt.deform = false;
    opt.scale = 'auto';
    opt.component = [];
    opt.colorbar = '';
    opt.shownodes = false;
    opt.elementlabels = true;
    opt.nodelabels = true;
    opt.new = true;
    opt.figure = [];
    opt.axes = []; 
    opt.plot = {};
    opt.patch = {};
    opt.hold = true;
    
    % When using data, disable the labels by default and do not
    % create a new figure with each call
    if ~isempty(opt.data);
        opt.elementlabels = false;
        opt.nodelabels = false;
        opt.new = false;
    end
    
    % Collect options supplied by the user
    opt = gatherUserOptions(opt, varargin{:}); 

    % Set/Create the figure handle
    if opt.new;
        figure(); hold on;   
    elseif ishandle(opt.figure)
        figure(opt.figure); hold on;
    end
    
    % Set the axes handle
    if ishandle(opt.axes);
        axes(opt.axes);
    end
    
    % Case when automatic scale factor for deformation is desired
    if ~isempty(opt.data) && opt.deform && strcmpi(opt.scale,'auto');
        n = max(max(floor(log10(abs(opt.data)))));
        opt.scale = 10^(-sign(n)*(abs(n)-1));  
    end
    
    % Account for text input for component option
    if ischar(opt.component);
        switch lower(opt.component);
            case 'x'; opt.component = 1;
            case 'y'; opt.component = 2;
            case 'z'; opt.component = 3;
            otherwise
                error('FEplot:parse_input','Un-reconginzed option, %s, for component option.', opt.component);
        end
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
            dof = elem.getDof();
            y = opt.data(dof);
        else
            y = zeros(size(x));
        end
        
        % Make sure y is a column vector
        if ~iscolumn(y); y = y'; end
        
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

function h = plot1D_vector(obj, opt)
    %PLOT1D_vector create a plot for 1D with multiple dofs per node

    % Initialize the x and y values
    h = zeros(obj.n_elements,1);
    
    % Generate the graph
    for e = 1:obj.n_elements;

        % The current element
        elem = obj.element(e);
        
        % The node positions
        x = elem.nodes(:,1);
        y = zeros(size(x));
        
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
                % Adjust the nodal values
                x = x + opt.scale*zz(:,1);
                y = z + opt.scale*zz(:,2);
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
            h(e) = patch(x, y, z,'EdgeColor','interp');
        end
    end
    
    % Create the colorbar        
    if ~isempty(opt.colorbar) && ischar(opt.colorbar);
        cbar = colorbar;
        set(get(cbar,'YLabel'),'String',opt.colorbar);
    end    

end

function h = plot2D_scalar(obj, opt)
    %PLOT2D_SCALAR create a 2D plot with scalar values

    % Initialize the x and y values
    h = zeros(obj.n_elements,1);
    
    % Generate the graph
    for e = 1:obj.n_elements;

        % The current element
        elem = obj.element(e);
        
        % Apply the plotting order
        if ~isempty(elem.node_plot_order);
            idx = elem.node_plot_order;
        else
            idx = 1:elem.n_dof;
        end

        % Create the graph
        if isempty(opt.data);
            h(e) = patch(elem.nodes(idx,1), elem.nodes(idx,2), 'w');
        else
            dof = elem.getDof();
            h(e) = patch(elem.nodes(idx,1), elem.nodes(idx,2),...
                opt.data(dof(idx)), 'EdgeColor','k');
        end
    end
end

function h = plot2D_vector(obj, opt)
    %PLOT2D_vector create a 2D plot for vector spaces
  
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
            dof = elem.getDof();
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
                if ~isempty(opt.component) && opt.component == 1 || strcmpi(opt.component, 'x');
                    x = x + opt.scale*zz(:,1);
                elseif ~isempty(opt.component) && opt.component == 2 || strcmpi(opt.component, 'y');
                    y = y + opt.scale*zz(:,2);
                else
                    x = x + opt.scale*zz(:,1);
                    y = y + opt.scale*zz(:,2);
                end
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
            h(e) = patch(x, y, z,'EdgeColor','interp');
        end
    end
end

function h = plot3D_scalar(obj, opt)
    %PLOT_SCALAR create a plot of scalar data

    % Initialize the x and y values
    h = zeros(obj.n_elements,1);
    
    % Generate the graph
    for e = 1:obj.n_elements;

        % The current element
        elem = obj.element(e);

        % Plot with data
        if ~isempty(opt.data)
            dof = elem.get_dof();
            h(e) = patch('Vertices',elem.nodes,'Faces',elem.side_dof,...
                'FaceVertexCData',opt.data(dof));
            
        % Plot without data    
        else
            h(e) = patch('Vertices',elem.nodes,'Faces',elem.side_dof);
        end
    end
    
    % Apply settings
    set(h,'FaceColor','none','EdgeColor','interp','Marker','.',...
        'MarkerFaceColor','flat');
    view(3);
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
    if ~isempty(opt.patch);
        set(h, opt.patch{:}); 
    end
    
    % Create the colorbar        
    if ~isempty(opt.colorbar) && ischar(opt.colorbar);
        cbar = colorbar;
        set(get(cbar,'YLabel'),'String',opt.colorbar);
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
    switch obj.options.type;
        
        % Continous finite element labels
        case 'CG'; 

            % Extract the unique nodes
            nodes = unique(obj.map.node,'rows','stable');
            
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
