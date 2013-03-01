function work

nx = 3;
ny = 3;
[nodes,node_map] = buildGrid(0,1,0,1,nx,ny);

spmd
    n_nodes = length(nodes);
    id = reshape(1:n_nodes,nx+1,ny+1)';
    n_elem = nx*ny;
    elements = cell(n_elem,1);
    k = 0;
    for j = 1:ny;
        for i = 1:nx;
            k = k + 1;
            elem_map(k,:) = reshape(id(i:i+1,j:j+1),4,1);
        end
    end
 
    codist = codistributor1d(1, codistributor1d.unsetPartition, size(elem_map));
    dist_elem_map = codistributed(elem_map, codist);
    
    part = [0,cumsum(codist.Partition)];
    e = ones(codist.Partition,4,'uint32');
    for i = 1:length(part)-1;
        e(part(i)+1:part(i+1),:) = i;
    end
    
    proc = reshape(e,numel(e),1);
    elem = reshape(elem_map, numel(elem_map),1);

%     local_elem_map = getLocalPart(elem_map);
%     local_elem_map_id = globalIndices(elem_map,1);

     node_codist = getCodistributor(nodes);
     
     % Nodes owned by this process
     local_node_id = globalIndices(nodes,1);
     

     [C,ia,ib] = intersect(local_node_id, elem,'stable');
     
     elem(ia)
     proc(ia)

     
     
     
     
     
    

    
    

end







    