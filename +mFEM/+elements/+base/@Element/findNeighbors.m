function findNeighbors(obj)
    %FINDNEIGHBORS Locates neighboring elements based on node coordinates
    %   This action is called by the Mesh class but is located within the
    %   element class for speed reasons.
    %
    % Syntax
    %   findNeighbors()
    %
    % Description
    %   findNeighbors() using the node coordintes neighboring elements are
    %   identified and the sides data structures are updated.
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
    
	% Loop over each of the elements
    for i = 1:length(obj);
        elem = obj(i);                           % current element
        elem_coord = elem.nodes.getCoord();      % current node coordinates
        neighbors = elem.nodes.getParents(elem); % neighbor elements

        % Loop through the element sides
        for s = 1:elem.n_sides;
            
            % Move to the next side if the neighbor is already defined
            if ~isempty(elem.sides(s).neighbor); 
                continue; 
            end
            
            % Get the coordinates of the current side
            elem_side = sort(elem_coord(:,elem.side_ids(s,:)),2);
            
            % Loop through the neighbor elements
            for j = 1:length(neighbors);
                neigh = neighbors(j);                   % neighbor element
                neigh_coord = neigh.nodes.getCoord();   % neighbor coords
                
                % Loop through neighbor sides
                for n = 1:neigh.n_sides;
                    neigh_side = sort(neigh_coord(:,neigh.side_ids(n,:)),2);

                    % If the sides are the same, update the structure
                    if isequal(elem_side, neigh_side);
                        elem.sides(s).neighbor = neigh;
                        elem.sides(s).neighbor_side = n;
                        elem.sides(s).on_boundary = false;
                        neigh.sides(n).neighbor = elem;
                        neigh.sides(n).neighbor_side = s;
                        neigh.sides(n).on_boundary = false;
                    end
                end
            end
        end

        % Update the element boundary flag
        idx = [elem.sides.on_boundary];
        elem.nodes(elem.side_ids(idx,:)).setBoundaryFlag();
        if ~all(idx); elem.on_boundary = true; end
    end
end