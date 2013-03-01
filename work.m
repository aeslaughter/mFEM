function work

opt.debug = true;
nx = 100;
ny = 50;
[nodes, node_map] = buildGrid(0,1,0,1,nx,ny);

spmd
    E.n_elem = nx*ny;
    E.elem_map = zeros(E.n_elem,4);
    
    n_nodes = length(nodes);
    id = reshape(1:n_nodes,nx+1,ny+1);

    k = 0;
    for j = 1:ny;
        for i = 1:nx;
            k = k + 1;
            E.elem_map(k,[1,2,4,3]) = reshape(id(i:i+1,j:j+1),4,1);
        end
    end
    
    E.codist = codistributor1d(1, codistributor1d.unsetPartition, size(E.elem_map));
    elem_map = codistributed(E.elem_map, E.codist);
%     local = getLocalPart(nodes);
%     node_id = globalIndices(nodes,1);
%     elem_id = globalIndices(elem_map,1);
    

    [info,ghost,ghost_id] = sendRecieveNodes(E,nodes);
    
    if opt.debug;
        disp(['Nodes ', mat2str(info.send.id), ' sent to ', mat2str(info.send.lab)]);
        disp(['Nodes ', mat2str(info.rec.id), ' are recieved from labs ', mat2str(info.rec.lab)]);
    end

end
end




function [info,ghost,ghost_id] = sendRecieveNodes(E, nodes)

    E.proc_map = buildElemProcMap(E);
    
    node_id = globalIndices(nodes,1);
    local = getLocalPart(nodes);

    send = struct('id',[], 'lab',[]);
    rec  = struct('id',[], 'lab',[]);

    for i = 1:E.n_elem;
       if labindex ~= E.proc_map(i);
           Lia = ismember(E.elem_map(i,:), node_id);
           send.id = [send.id, E.elem_map(i,Lia)];
           send.lab = [send.lab,repmat(E.proc_map(i),1,sum(Lia))];
           
       elseif labindex == E.proc_map(i);
           Lia = ~ismember(E.elem_map(i,:), node_id);
           rec.id = [rec.id, E.elem_map(i,Lia)];
           rec.lab = [rec.lab,repmat(E.proc_map(i),1,sum(Lia))];
       end
    end

    [send.id,ia] = unique(send.id);
    send.lab = send.lab(ia);

    [rec.id,ia] = unique(rec.id);
    rec.lab = rec.lab(ia);
    
    info = struct('send', send, 'rec', rec);


    [~,idx] = intersect(node_id, send.id, 'stable');
    for i = 1:length(idx)
        labSend(local(idx(i)), send.lab(i), 1);
        labSend(send.id(i), send.lab(i), 2);
    end
    
    N = length(rec.id);
    ghost = cell(N,1);
    ghost_id = zeros(N,1);
    for i = 1:N
        ghost{i} = labReceive;%(rec.lab(i));
        ghost_id(i,1) = labReceive;
    end
end

function proc_map = buildElemProcMap(E)
    % Build elem_proc_map
    proc_map = zeros(E.n_elem,1);
    part = [0,cumsum(E.codist.Partition)];
    for i = 1:length(part)-1;
        proc_map(part(i)+1:part(i+1)) = i;
    end
end





    