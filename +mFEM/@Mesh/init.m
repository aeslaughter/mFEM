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

    % Initilize boundary map
%     no_boundary_ids = getpref('MFEM_PREF','NO_BOUNDARY_IDS');
%     obj.boundary_tag
%     spmd
%         codist = getCodistributor(node_map);
%         gsize = [n,no_boundary_ids];
%         codist = codistributor(1,codist.Partition,gsize);
%         boundary_map = codistributed.false(gsize,codist);
%         
%         local = getLocalPart(boundary_map);
%     end
%     
    
end