function plot(obj, data, varargin)

    opt.data = data;
    
    if isempty(opt.data);
        opt.labelelements = true;
        opt.labelnodes = true;   
    else
        opt.labelelements = false;
        opt.labelnodes = false;
    end
    opt = gatherUserOptions(opt,varargin{:});

    p = zeros(length(obj.elements),1);

    figure; hold on;
    elements = gather(obj.elements);
    
    for i = 1:length(elements);
        local = elements{i};
        for e = 1:length(local);
            elem = local(e);
            vert = elem.nodes.getCoord();
            face = elem.side_ids;
            p(i) = patch('Vertices',vert,'Faces',face);

            if opt.labelelements;
                cntr = num2cell(mean(vert,1));
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
    end
end