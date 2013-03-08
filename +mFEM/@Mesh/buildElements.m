function [elements,nodes] = buildElements(obj,nodes) 
    %BUILDELEMENTS (protected) Create elements from element map and nodes
    %
    % Syntax
    %   [elements,nodes] = buildElements(nodes)
    %
    % Description
    %
    %
    
    % Get variables for use
    elem_map = obj.elem_map;
    node_map = obj.node_map;
    
    % Get the element type (todo) this will allow for mixed meshes
    elem_type = obj.elem_type;
    type = gather(elem_type(1));type = type{1};
    
    % Serial case
    if matlabpool('size') == 0;
        spmd
            n_elem = size(elem_map,1); % no. of elements
            elements(n_elem) = feval(['mFEM.elements.',type]); 
            elements.init(1:n_elem, nodes(gather(elem_map)));
        end
        return % done
    end

    % Parallel case
    spmd
        % Extract local part of the element map 
        e_map = getLocalPart(elem_map);
        e_id  = globalIndices(elem_map, 1); % global element no.

        % Extract all the nodes that are needed on this processor
        [no, no_id] = getOffLabNodes(e_map, node_map, nodes);
        n = size(e_map,1);
        [~,Locb] = ismember(reshape(e_map,1,numel(e_map)),no_id);
        e_map = reshape(Locb, size(e_map));
        
        % Create the elements
        type = elem_type{1}; %(todo)
        elements(n,1) = feval(['mFEM.elements.',type]);
        elements.init(e_id, no(e_map));  
    end
end

function [local,local_id] = getOffLabNodes(e_map, node_map, nodes)
    %GETOFFLABNODES gets nodes that are needed from other labs
    
    % Extract partition information from the nodes array
    n_dist = getCodistributor(node_map);
    part = cumsum(n_dist.Partition);
    local = nodes;
    local_id = globalIndices(node_map,1);  % local node ids
    limit = getpref('MFEM_PREF','LABSEND_LIMIT');
   
    % Build a map that indicates the lab location for each node
    no = unique(e_map);         % all nodes needed by this lab
    loc = ones(size(no));       % initilize location map
    for i = 2:numlabs
        loc(no > part(i-1)  & no <= part(i)) = i;
    end
    
    % Seperate the nodes needed that are not on the current lab
    ix = loc~=labindex;    % excludes nodes already on this lab
    node_ids = no(ix);      % off lab node ids needed
    lab_map = loc(ix);      % locations of the nodes needed
    lab = unique(loc(ix));  % list of labs that will need to be accessed
  
    % Loop through each of the labs containing nodes that are needed
    for i = 1:length(lab);
        % Send a request to the lab for nodes based on the id
        labSend(node_ids(lab_map==lab(i)),lab(i));
        %disp(['Sending request to lab ', num2str(lab(i)), ' for nodes ', mat2str(node_ids(lab_map==lab(i)))]);
    end
    labBarrier; % finish all sends before continuing
     
    % Recieve the request for nodes
    k = 1;
    nid = {};
    while labProbe % continue to receive requests until there is none
        [nid{k},proc(k)] = labReceive;    
        %disp(['Lab ', num2str(proc(k)), ' requested nodes ', mat2str(nid{k})]);
        k = k + 1;    
    end
    labBarrier; % finish recieving requests for nodes before continuing
    
    % Send a package (struct) containing the nodes and global ids
    for k = 1:length(nid);
        [~,idx] = intersect(local_id,nid{k});

        s = [0:limit:length(idx),length(idx)];
        for i = 1:length(s)-1;
            ix = idx(s(i)+1:s(i+1)); 
            pkg.id = local_id(ix);
            pkg.nodes = local(ix);
            labSend(pkg,proc(k));
            %disp(['Lab ', num2str(labindex), ' sent nodes ', mat2str(pkg.id), ' to lab ', num2str(proc(k))]);
        end
    end
    labBarrier; % finish sending nodes before continuing
     
    % Recieve the package and append local nodes
    while labProbe
        [pkg,p] = labReceive;           % recieve the package
        local = [local; pkg.nodes];     % append recieved nodes
        local_id = [local_id, pkg.id];  % append recieved node ids
        %disp(['Lab ', num2str(labindex), ' recieved nodes ', mat2str(pkg.id), ' from lab ', num2str(p)]);
    end
    labBarrier;
end