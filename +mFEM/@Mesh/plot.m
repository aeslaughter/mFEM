function plot(obj, data, varargin)

    % Check if system is initialized
%     if ~obj.initialized;
%         error('Mesh:plot:NonInitializedMesh',...
%             'The Mesh object must be initialized for plotting');
%     end
 
    opt.data = data;
    opt.shownodes = false;
    
    if isempty(opt.data);
        opt.labelelements = true;
        opt.labelnodes = true;   
    else
        opt.labelelements = false;
        opt.labelnodes = false;
    end
    opt = gatherUserOptions(opt,varargin{:});


    elements = gather(obj.elements);

    hfig = gcf;
    set(hfig,'NextPlot','add');
    k = 0;
    for i = 1:length(elements);
        local = elements{i};
        for e = 1:length(local);
            k = k + 1;
            elem = local(e);
            vert = elem.nodes.getCoord();

            face = elem.side_ids;
            p(k) = patch('Vertices',vert','Faces',face);
            
            if ~isempty(opt.data);
                opt.data
                set(p(k),'FaceVertexCdata',opt.data);
            end

            if opt.labelelements;
                cntr = num2cell(mean(vert',1));
                text(cntr{:}, num2str(elem.id),'FontSize',14,...
                'BackgroundColor','k','FontWeight','Bold',...
                'Color','w','HorizontalAlignment','center');
            end

            if opt.labelnodes;
                for j = 1:length(elem.nodes);
                    n = elem.nodes(j);
                    loc = num2cell(n.coord);
                    text(loc{:},num2str(n.id),'FontSize',10,'Color','w',...
                        'BackgroundColor','b','HorizontalAlignment','center');
                end
            end
        end
        
        if opt.shownodes;
            set(p,'Marker','o','MarkerEdgeColor','k');
            
        end
    end
end