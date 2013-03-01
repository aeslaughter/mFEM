function work

nx = 3;
ny = 3;
[nodes,node_map] = buildGrid(0,1,0,1,nx,ny);

spmd
    n_nodes = length(nodes);
    id = reshape(1:n_nodes,nx+1,ny+1);
    n_elem = nx*ny;
    elements = cell(n_elem,1);
    k = 0;
    for j = 1:ny;
        for i = 1:nx;
            k = k + 1;
            elem_map(k,[1,2,4,3]) = reshape(id(i:i+1,j:j+1),4,1);
        end
    end
    
    codist = codistributor1d(1, codistributor1d.unsetPartition, size(elem_map));
    dist_elem_map = codistributed(elem_map, codist);

    [no_send, no_send_lab] = locateSendNodes(globalIndices(nodes,1), elem_map, codist);
    disp(['Nodes ', mat2str(no_send), ' are needed by labs ', mat2str(no_send_lab)]);
    [no_rec, no_rec_lab] = locateRecieveNodes();
    
    
    
    % Build elem_proc_map
%     elem_proc_map = zeros(size(elem_map));
%     part = [0,cumsum(codist.Partition)];
%     for i = 1:length(part)-1;
%         elem_proc_map(part(i)+1:part(i+1),:) = i;
%     end
% 
%     elem_proc_map = reshape(elem_proc_map, numel(elem_proc_map), 1);
%     elem_linear_map = reshape(elem_map, numel(elem_map), 1)
    
    
%     local_node_id = globalIndices(nodes,1);
%     local_elem_id = globalIndices(dist_elem_map,1);

    % Get nodes needed by other processors
%     idx = true(n_elem,1);
%     idx(local_elem_id) = false;
%     off_proc_elem_nodes = unique(reshape(elem_map(idx,:),sum(idx)*4,1));
    
    [~,~,ib] = intersect(local_node_id, elem_linear_map, 'stable');
    
    no = elem_linear_map(ib);
    loc = elem_proc_map(ib)~=labindex;
    no = no(loc);
    
    disp(['Nodes ', mat2str(no), ' are needed by labs ', mat2str(elem_proc_map(loc))]);

    
    


    
%     local_elem_map = getLocalPart(dist_elem_map);
%     no = unique(reshape(local_elem_map,numel(local_elem_map),1));
%     
%     no = no(no < local_node_id(1) | no > local_node_id(end));
    
    % Locate the nodes needed


end

function [no,lab] = locateSendNodes(node_id, elem_map, codist)

    % Build elem_proc_map
    n_elem = size(elem_map,1);
    elem_proc_map = zeros(n_elem,1);
    part = [0,cumsum(codist.Partition)];
    for i = 1:length(part)-1;
        elem_proc_map(part(i)+1:part(i+1)) = i;
    end
    
    no = [];
    lab = [];
    for i = 1:n_elem;
       if labindex ~= elem_proc_map(i);
            [C,ia,ib] = intersect(node_id, elem_map(i,:));
            no = [no,C];
            lab = [lab,repmat(elem_proc_map(i),1,length(C))];
       end
    end
    
[no,ia] = unique(no);
lab = lab(ia);

%     elem_proc_map = reshape(elem_proc_map, numel(elem_proc_map), 1);
%     elem_linear_map = reshape(elem_map, numel(elem_map), 1)
%     
%     
%     local_node_id = globalIndices(nodes,1);
% %     local_elem_id = globalIndices(dist_elem_map,1);
% 
%     % Get nodes needed by other processors
% %     idx = true(n_elem,1);
% %     idx(local_elem_id) = false;
% %     off_proc_elem_nodes = unique(reshape(elem_map(idx,:),sum(idx)*4,1));
%     
%     [~,~,ib] = intersect(local_node_id, elem_linear_map, 'stable')
%     
%     no = elem_linear_map(ib);
%     loc = elem_proc_map(ib)~=labindex;
%     no = no(loc);
%     
%     disp(['Nodes ', mat2str(no), ' are needed by labs ', mat2str(elem_proc_map(loc))]);






    