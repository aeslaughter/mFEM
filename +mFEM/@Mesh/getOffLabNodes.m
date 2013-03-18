function local = getOffLabNodes(elem_map, node_map, nodes)
    %GETOFFLABNODES gets nodes that are needed from other labs
    %   Each element contains the nodes that are associated with the
    %   element, if the nodes are stored on other processors they must be
    %   retrieved. This function sends requests to labs for nodes, which in
    %   turn send back the nodes (copies).
    %
    %   This function must be run inside an spmd block.
    %
    % Syntax
    %   no = getOffLabNodes(e_map, node_map, nodes)
    %
    % Description
    %   no = getOffLabNodes(e_map, node_map, nodes) gets
    %   the necessary nodes: e_map is the codistributed elem_map,
    %   node_map is the codistributed nodal map, and nodes is the
    %   Composite vector of node objects. It returns an array of node
    %   objects ready for insertion into the Element init() method.
    %
    %----------------------------------------------------------------------
    %  mFEM: A Parallel, Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2013 Andrew E Slaughter
    % 
    %  This program is free software: you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation, either version 3 of the License, or
    %  (at your option) any later version.
    % 
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    % 
    %  You should have received a copy of the GNU General Public License
    %  along with this program. If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
    
    % Extract partition information from the nodes array
    n_dist = getCodistributor(node_map);
    e_map = getLocalPart(elem_map);
    part = cumsum(n_dist.Partition);
    local = nodes;
    local_id = [nodes.id]; % local node ids
    limit = getpref('MFEM_PREF','LABSEND_LIMIT');
   
    % Build a map that indicates the lab location for each node
    no = unique(e_map);         % all nodes needed by this lab
    loc = ones(size(no));       % initilize location map
    for i = 2:numlabs
        loc(no > part(i-1)  & no <= part(i)) = i;
    end
    
    % Seperate the nodes needed that are not on the current lab
    ix = loc~=labindex;     % excludes nodes already on this lab
    node_ids = no(ix);      % off lab node ids needed
    lab_map = loc(ix);      % locations of the nodes needed
    lab = unique(loc(ix));  % list of labs that will need to be accessed
  
%     disp(sprintf('Data is available: %d',labProbe));
    
    % Loop through each of the labs containing nodes that are needed
    for i = 1:length(lab);
        % Send a request to the lab for nodes based on the id
        labSend(node_ids(lab_map==lab(i)),lab(i),100);
        disp(['Lab ', num2str(labindex), ' sent request to lab ', num2str(lab(i)), ' for nodes: ', mat2str(node_ids(lab_map==lab(i)))]);
    end
    labBarrier;
    
    % If a request is found, send out those nodes
    while labProbe('any',100) % continue to receive requests until there is none
        [nid,proc] = labReceive('any',100);  
        [~,idx] = intersect(local_id,nid);
        
        % Limite the package size
        s = [0:limit:length(idx),length(idx)];
        for i = 1:length(s)-1;
            ix = idx(s(i)+1:s(i+1));
            labSend(local(ix),proc,200);
        end
        
        disp(['Lab ', num2str(proc), ' requested nodes from lab ',num2str(labindex),': ', mat2str(nid)]);
        disp(['Lab ', num2str(labindex), ' sent nodes to lab ', num2str(proc),': ', mat2str(nid)]);
    end
    labBarrier;
     
    % Recieve the package and append local nodes
%     
    for i = 1:length(lab);
        while labProbe(lab(i),200);
            [ghst,p] = labReceive(lab(i),200);  % recieve the package
            local = [local; ghst];     % append recieved nodes
        end
        disp(['Lab ', num2str(labindex), ' recieved nodes from lab ', num2str(p)]);
    end
    labBarrier;
    
    % Convert global node connectivity element map to local indices 
    [~,Locb] = ismember(reshape(e_map,1,numel(e_map)),[local.id]);
    e_map = reshape(Locb, size(e_map));

    % Ensure that nodes are organized correctly (single element)
    local = local(e_map);
    if ~isequal(size(e_map),size(local));
        local = local';
    end
end