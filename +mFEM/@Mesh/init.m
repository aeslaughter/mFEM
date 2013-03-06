function init(obj)

    %
    nodes    = obj.nodes;
    elements = obj.elements;
    node_map = obj.node_map;
    n_nodes  = size(obj.node_map,1);
    obj.n_nodes = n_nodes;
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
    
    space = obj.options.space;
    spmd
        nodes.setDof(space);
    end

    obj.initialized = true;
end