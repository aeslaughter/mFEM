function update(obj)
    %UPDATE Updates nodes of elements
    %   The nodes stored in the elements are copies of the nodes stored in
    %   within the Mesh:nodes property. This function must be called after
    %   any changes to the nodes are made (i.e., after addBoundary or
    %   addSubdomain). This is only a problem for elements that have nodes
    %   from other processors.
    %
    % Syntax
    %   update()
    %
    % Description
    %   update() updates the nodes on the elements in parallel operation
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
     
    % Perform parallel node update
    elem_map = obj.elem_map;
    node_map = obj.node_map;
    elements = obj.elements;
    nodes = obj.nodes;   

    % No communication is required in serial
    if matlabpool('size') == 0;
        spmd
            elements.update();
        end
        return
    end
    
    % Parallel requires off processor nodes to be copied
    spmd
        % Extract all the nodes that are needed on this processor
        no = obj.getOffLabNodes(elem_map, node_map, nodes);
        
        % Update the elements nodes
        [m,n] = size(no);
        no = mat2cell(no,ones(m,1),n);
        [elements.nodes] = no{:};
        labBarrier;
        
        % Update/apply the element tags
        elements.update();
    end
end
    