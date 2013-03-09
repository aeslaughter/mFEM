function init(obj)
    %INIT (protected) Initializes the mFEM.Mesh object for use
    %   The setup process does three main tasks: (1) locating
    %   neighboring elements, (2) sets the degrees-of-freedom for nodes,
    %   and (3) sets up the various tag maps for identifing boundaries and
    %   subdomains.
    %
    % Syntax
    %   init()
    %
    % Description
    %   init() prepares the mFEM.Mesh object for use. by locating the
    %   neighboring elements, setting the degrees-of-freedom for each node
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

    % Build the nodes
    nodes = obj.buildNodes();
    
    % Build the elements
    [elements,nodes] = obj.buildElements(nodes);
    
    % Set no. of nodes and elements
    obj.n_nodes = size(obj.node_map,1);
    obj.n_elements = size(obj.elem_map,1);
    
    % Display message time for finding neighbors, if desired
    if obj.options.time;
        ticID = tMessage('Locating neighbor elements...');
    end

    % Perform neighbor location, in parallel. This is done within the
    % element class for speed reasons.
    spmd
        elements.findNeighbors();
    end

    % Complete message time message
    if obj.options.time;
        tMessage(ticID);
    end
    
    % Set the nodal degrees-of-freedom, in parallel
    space = obj.options.space;
    n_nodes = obj.n_nodes;
    spmd
       ndof = buildCodistributedPartition(sum([nodes.n_dof]));
       strt = [0,cumsum(ndof)] + 1;  
       nodes.setDof(strt(labindex));
    end
    
    % Set the total degrees-of-freedom for the mesh
    obj.n_dof = sum(ndof{1});
    
    % Update the Mesh object variables
    obj.nodes = nodes;
    obj.elements = elements;    
    
    % The following code initializes the tag maps for the element and 
    % nodes, in parallel. Each map is a codistributed sparse logical 
    % matrix.
    
    % Extract the maximum number of boundary ids from the preferences
    n_col = getpref('MFEM_PREF','MAX_NUM_TAGS');
    no_dist = obj.node_map_codist;  % codistribitor of nodes (need Partition)
    no_gsize = [obj.n_nodes, n_col];% global size of node tag map
    el_dist = obj.elem_map_codist;  % codistribitor of elements (need Partition)
    el_gsize = [obj.n_elements, n_col]; % global size of element tag map
    
    % Begin parallel operations
    spmd 
        % Create the sparse logical matrix for node tags
        no_dist = codistributor1d(1, no_dist.Partition, no_gsize);
        node_tag_map = logical(sparse(no_gsize(1),no_gsize(2),no_dist));
        
        % Create the sparse logical matrix for element tags
        el_dist = codistributor1d(1, el_dist.Partition, el_gsize);
        elem_tag_map = logical(sparse(el_gsize(1),el_gsize(2),el_dist));
    end
    
    % Update the mesh object with initalized tag maps
    obj.node_tag_map = node_tag_map;
    obj.elem_tag_map = elem_tag_map;
end