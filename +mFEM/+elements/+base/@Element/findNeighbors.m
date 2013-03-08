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
%         m = elem.side_ids(idx,:);
%         in = num2cell(true(numel(m),1));
%         elem.nodes(m).on_boundary
%         [elem.nodes(m).on_boundary] = in{:};
        elem.nodes(elem.side_ids(idx,:)).setBoundaryFlag();
        if ~all(idx); elem.on_boundary = true; end
    end
% ISMEMBER METHOD: very slow (half as above)          
%             for i = 1:length(obj);
%                 elem = obj(i);
%                 neighbors = elem.nodes.getParents(elem);
%                 n = size(elem.side_ids,2);
%                 
%                 for j = 1:length(neighbors);
%                      neigh = neighbors(j);
%                      
%                      nse = neigh.getSideElements();
%                      if any(nse == elem); continue; end;
%                     idx = ismember(elem.smap.map,neigh.smap.map,'rows');
%                     h = histc(elem.smap.id(idx),elem.smap.uid);
%                     s1 = elem.smap.uid(h >= n);
%                     if isempty(s1); continue; end
%                     
%                     idx = ismember(neigh.smap.map,elem.smap.map,'rows');
%                     h = histc(neigh.smap.id(idx),neigh.smap.uid);
%                     s2 = neigh.smap.uid(h >= n);
% 
%                     elem.sides(s1).neighbor = neigh;
%                     elem.sides(s1).neighbor_side = s2;
%                     elem.sides(s1).on_boundary = false;
%                     neigh.sides(s2).neighbor = elem;
%                     neigh.sides(s2).neighbor_side = s1;  
%                     neigh.sides(s2).on_boundary = false;
%                 end
%             end
end