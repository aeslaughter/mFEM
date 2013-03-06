function init(obj)

    nodes = obj.nodes;
    space = obj.options.space;
    elements = obj.elements;
    obj.n_nodes = size(obj.node_map,1);
    obj.n_elements = size(obj.elem_map,1);
    
    % Display wait message
    if obj.options.time;
        ticID = tMessage('Locating neighbor elements...');
    end

    spmd
        elements.findNeighbors();
    end

    % Complete message
    if obj.options.time;
        tMessage(ticID);
    end;
    
    spmd
    %     codist = getCodistributor(node_map);
    %     codist = codistributor1d(1,codist.Partition,[obj.n_nodes,1]);  
    %     boundary_map = codistributed.false(obj.n_nodes,1,codist);
        nodes.setDof(space);
    end
end