function findNeighbors(obj)
% BRUTE FORCE: slow (fastest I can acheive)     
    for i = 1:length(obj);
        elem = obj(i);
        elem_coord = elem.nodes.getCoord();
        neighbors = elem.nodes.getParents(elem);

        for s = 1:elem.n_sides;
            if ~isempty(elem.sides(s).neighbor); 
                continue; 
            end
            elem_side = sort(elem_coord(:,elem.side_ids(s,:)),2);
            for j = 1:length(neighbors);
                neigh = neighbors(j);
                neigh_coord = neigh.nodes.getCoord();
                for n = 1:neigh.n_sides;
                    neigh_side = sort(neigh_coord(:,neigh.side_ids(n,:)),2);

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
        
        idx = [elem.sides.on_boundary];
        elem.nodes(elem.side_ids(idx,:)).setBoundaryFlag();
        if ~all(idx); elem.on_boundary = true; end
    end
    
    % Change parents from handles to ids
%     for i = 1:length(obj);
%         obj(i).nodes.resetParents();
%     end
end