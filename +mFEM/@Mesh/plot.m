function h = plot(obj, data, varargin)
    %PLOT Creates graphs for finite element data.
    %
    % Syntax
    %  plot([])
    %  plot(data)
    %  plot(...,'PropertyName',<PropertyValue>)
    %
    % Descrtiption
    %   Fplot([]) creates a plot of the mesh object (obj) and includes that
    %   includes the element and degree of freedom labels.
    %
    %   plot(data) creates a plot for the mesh object (obj) and the data 
    %   supplied, the data must be a column vector of with the number of values 
    %   equal to the number of degrees of freedom for the mesh.
    %
    %   plot(...,'PropertyName',PropertyValue,...) allows for the
    %   plot created to be customised, for a list of properties below. Note 
    %   that all true/false may be input using a flag instead of the pairing. 
    %   For example the following two calls to FEplot are the same:
    %       FEplot(obj,data,'-Deform');
    %       FEplot(obj,data,'Deform',true);
    %   This flag style input simply reverses the state from the default.
    %
    % plot Property Descriptions
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
    %       true | {false}
    %       Toggles the display of plotted dots at the nodes
    %
    %  ShowElements
    %       true | {false}
    %       Toggles the display of the element grid lines
    %
    %  ElementLabels
    %       true (default w/o data) | false (default w/ data)
    %       Toggles the appearence of element number labels.
    %       
    %  NodeLabels
    %       true (default w/o data) | false (default w/ data)
    %       Toggles the appearence of node number labels.
    %
    %  New
    %       true (default w/o data) | false (default w/ data)
    %       Toggles the creation of a new figure, if it is set to false then
    %       gcf is used, unless the handle is specified via FigureHandle
    %
    %  Axes
    %       axes handle
    %       Add the plot to the user specified axes handle. 
    %
    %  Figure
    %       figure handle
    %       Add the plot to the user specified figure handle. Note, NewFigure
    %       takes presidence over FigureHandle.
    %
    %  Patch
    %       cell array
    %       Use this to pass commands directly to patch (2D & 3D), for example
    %           FEplot(...,'Patch',{'FaceColor','k'});
    %       executes the following:
    %           patch(...,'FaceColor','k');
    %       This is the last command to be applied to the patch, thus it will
    %       override other patch related settings.
    %
    %   Slice
    %       numeric matrix

    % See also Mesh
    %
    %----------------------------------------------------------------------
    %  mFEM: A Parallel, Object-Oriented MATLAB Finite Element Library
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
    %  along with this program. If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------

    % Collect the input 
    opt = parseInput(data, varargin{:});

%     elements = obj.getElements();
    
%     for i = 1:length(elements);
%         
%     
% %     order = elements.getNodePlotOrder();   
%     
    tmp = gather(obj.elem_map);
    for i = 1:length(tmp);
        faces(i,:) = tmp(i,[1,2,3,4]);
    end
    
    FV.faces = faces;
    FV.vertices = gather(obj.node_map);

    if size(FV.vertices,2) == 1 && isempty(opt.data);
        FV.vertices(:,2) = 0;
    elseif size(FV.vertices,2) == 1;
        FV.vertices(:,2) = opt.data;
    end

    space = obj.options.space;
    if (ischar(space) && strcmpi(space,'scalar')) || ...
       (isnumeric(space) && space == 1)
        FV.FaceVertexCData = opt.data;
        
    elseif (ischar(space) && strcmpi(space,'vector')) || ...
           (isnumeric(space) && space > 1);
       
        n = obj.n_dof/obj.n_nodes;
        opt.data = reshape(opt.data,n,numel(opt.data)/n)';
        
        if isempty(opt.component);
            data = mean(opt.data,2);
        else
            data = opt.data(:,opt.component);
        end
        FV.FaceVertexCData = data;
    end
    
    if isfield(FV,'FaceVertexCData') &&  opt.deform;
        if ~isempty(opt.component);
            c = opt.component;
            FV.vertices(:,c) = FV.vertices(:,c) + opt.scale*opt.data(:,c); 
        else
            FV.vertices = FV.vertices + opt.scale*opt.data;
        end
    end
    
    if isempty(FV.FaceVertexCData);
        h = patch(FV);
    else
        h = patch(FV,'FaceColor','interp');
    end
    
    applyPlotOptions(obj, h, opt);
%     
%     % Add the element labels
%     addElementLabels(obj, opt);
%        
%     % Add the nodel lables
%     addNodeLabels(obj, opt)
%     
%     % Resize the figure
%     box on;
end 

function opt = parseInput(data, varargin)
    %PARSEINPUT Function for parsing input data to plot function
    %
    % See FEPLOT

    % Define user properties
    opt.data = data; 
    opt.deform = false;
    opt.scale = 'auto';
    opt.component = [];
    opt.colorbar = '';
    opt.shownodes = false;
    opt.showelements = false;
    opt.elementlabels = true;
    opt.nodelabels = true;
    opt.new = true;
    opt.figure = [];
    opt.axes = []; 
    opt.plot = {};
    opt.patch = {};
    opt.hold = true;
    opt.slice = {};
    
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
    opt.component = convertComponent(opt.component);
end

function applyPlotOptions(obj,h,opt)
    % APPLY_PLOT_OPTIONS
    
    % Show the nodes as empty circles
    if opt.shownodes
        set(h, 'Marker','o', 'MarkerSize', 8, 'MarkerEdgeColor', 'auto',...
            'MarkerFaceColor', 'auto');
    end
    
    if opt.showelements
        set(h,'EdgeColor','interp');
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

function addElementLabels(obj, opt)
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

function addNodeLabels(obj, opt)
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
