function nodes = buildNodes(node_map)

    if matlabpool('size') == 0;
        node_map = gather(node_map);
        n = size(node_map,1);
        nodes(n) = mFEM.elements.base.Node();
        nodes.init(1:n,node_map);
        nodes = num2cell(nodes);
        return;
    end

% new method: 1000x1000 quad4 = 8 s
    spmd  
        local_map = getLocalPart(node_map);
        id = globalIndices(node_map,1);
        n = size(local_map,1);
        local(n,1) = mFEM.elements.base.Node();
        local.init(id,local_map);
        codist = getCodistributor(node_map);
        n_nodes = size(node_map,1);
        codist = codistributor1d(1, codist.Partition, [n_nodes,1]);
        nodes = codistributed.build(num2cell(local), codist);
    end
    
% old method: 1000x1000 quad4 = 29.2
%     spmd  
%         local_map = getLocalPart(node_map);
%         id = globalIndices(node_map,1);
%         n = size(local_map,1);
%         local = cell(n,1);
%         for i = 1:n;
%             local{i} = mFEM.elements.base.Node(id(i), local_map(i,:));
%         end
%         codist = getCodistributor(node_map);
%         n_nodes = size(node_map,1);
%         codist = codistributor1d(1, codist.Partition , [n_nodes,1]);
%         nodes = codistributed.build(local, codist);
%     end