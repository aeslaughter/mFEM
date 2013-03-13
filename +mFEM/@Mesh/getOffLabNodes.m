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
    local_id = globalIndices(node_map,1);  % local node ids
    limit = getpref('MFEM_PREF','LABSEND_LIMIT');
%     limit = 10;
   
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
%         disp(['Sending request to lab ', num2str(lab(i)), ' for nodes ', mat2str(node_ids(lab_map==lab(i)))]);
    end
    labBarrier; % finish all sends before continuing
     
    % Recieve the request for nodes
    k = 1;
    nid = {};
    while labProbe % continue to receive requests until there is none
        [nid{k},proc(k)] = labReceive;    
%         disp(['Lab ', num2str(proc(k)), ' requested nodes ', mat2str(nid{k})]);
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
    
    % Convert global node connectivity element map to local indices 
    [~,Locb] = ismember(reshape(e_map,1,numel(e_map)),local_id);
    e_map = reshape(Locb, size(e_map));

    % Ensure that nodes are organized correctly (single element)
    local = local(e_map);
    if ~isequal(size(e_map),size(local));
        local = local';
    end
end