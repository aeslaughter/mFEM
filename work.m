function work

x0 = 0; x1 = 1;
y0 = 0; y1 = 1;
xn = 10; yn = 10;

% tic;
% spmd
%     x = 1:(xn+1)*(yn+1);
%     id = reshape(x,xn+1,yn+1);
%     elem_map = zeros(yn*xn,4);
%     k = 0;
%     for j = 1:yn;
%         for i = 1:xn;
%             k = k + 1;
%             elem_map(k,[1,2,4,3]) = reshape(id(i:i+1,j:j+1),4,1);
%         end
%     end
%     codist = codistributor1d(1, codistributor1d.unsetPartition , size(elem_map));
%     elem_map = codistributed(elem_map, codist);
% end
% toc;


% tic;
% spmd
%     e = uint32(1:xn*yn);
%     elem_map(:,1) = e + ceil(e/(yn))-1;
%     elem_map(:,2) = elem_map(:,1) + 1;
%     elem_map(:,3) = elem_map(:,1) + xn + 2;
%     elem_map(:,4) = elem_map(:,1) + xn + 1;
%     codist = codistributor1d(1, codistributor1d.unsetPartition , size(elem_map));
%     elem_map = codistributed(elem_map, codist);
% end
% toc;
% elem_map

tic;
spmd
    codist = codistributor1d(1, codistributor1d.unsetPartition , [xn*yn,4]);
    part = [0, cumsum(codist.Partition)];
    e = uint32(part(labindex)+1:part(labindex+1));
    map(:,1) = e + ceil(e/(yn))-1;
    map(:,2) = map(:,1) + 1;
    map(:,3) = map(:,1) + xn + 2;
    map(:,4) = map(:,1) + xn + 1;
    elem_map = codistributed.build(map,codist);
end
toc;

tic;
spmd
    x = codistributed(x0:(x1-x0)/xn:x1);
    y = codistributed(y0:(y1-y0)/yn:y1);
    [X,Y] = ndgrid(x,y); 
    node_map = [reshape(X,numel(X),1), reshape(Y,numel(Y),1)]; 
end
toc;

tic;
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
toc;

tic;
spmd
    local_map = getLocalPart(elem_map);
    id = globalIndices(elem_map,1);
    n = size(local_map);
    
    nlocal= getLocalPart(nodes);
    nid = globalIndices(nodes,1);
    
    no = unique(local_map);
    
    ndist = getCodistributor(node_map);
    part = cumsum(ndist.Partition);
    
    loc = ones(size(no));
    for i = 2:numlabs
        loc(no>part(i-1)&no<=part(i)) = i;
    end
    
    ix = loc~=labindex;
    no = no(ix);
    loc = loc(ix);

    u = unique(loc);
    for i = 1:length(u);
        labSend(no(loc==u(i)),u(i));
        disp(['Sending request to lab ', num2str(u(i)), ' for nodes ', mat2str(no(loc==u(i)))]);
    end
    labBarrier;
    
    idx = {};
    k = 1;
    while labProbe
        [idx{k},L(k)] = labReceive;        
        disp(['Lab ', num2str(L(k)), ' requested nodes ', mat2str(idx{k})]);
        k = k + 1;
    end
    
    labBarrier;
    for k = 1:length(idx);
        [~,ib] = intersect(nid,idx{k});
        pkg.id = nid(ib);
        pkg.nodes = nlocal(ib);
        labSend(pkg,L(k));
        disp(['Lab ', num2str(labindex), ' sent nodes ', mat2str(ib), ' to lab ', num2str(L(k))]);
    end
    
    
    labBarrier;
    while labProbe
        [pkg,L(k)] = labReceive;
        pkg
%         [idx{k}] = labReceive;
%         disp(['Lab ', num2str(labindex), ' recieved nodes ', mat2str(idx{k}), ' from lab ', num2str(L(k))]);
        k = k + 1;
    end
    
end
toc;











% tic;
% spmd
%     elem_id = codistributed(uint32(0:nx*ny-1));
% 
%     e = getLocalPart(elem_id);
%     
% %     elem_map = zeros(length(e),4,'uint32');
%     map(:,1) = e + floor(e/ny) + 1;
%     map(:,2) = map(:,1) + 1;
%     map(:,3) = map(:,1) + 6;
%     map(:,4) = map(:,1) + 5;
% %     
% %     
% %     codist = getCodistributor(elem_id);
% %     part = codist.Partition;
%     codist = codistributor1d(1, codistributor1d.unsetPartition , [nx*ny,4]);
%     elem_map = codistributed.build(map,codist);
% end
% toc;


return;


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





    