function nodes = buildNodes(node_map)

    spmd  
        local_map = getLocalPart(node_map);
        id = globalIndices(node_map,1);
        n = size(local_map,1);
        local = cell(n,1);
        for i = 1:n;
            local{i} = mFEM.elements.base.Node(id(i), local_map(i,:));
        end
        codist = getCodistributor(node_map);
        n_nodes = size(node_map,1);
        codist = codistributor1d(1, codist.Partition , [n_nodes,1]);
        nodes = codistributed.build(local, codist);
    end
