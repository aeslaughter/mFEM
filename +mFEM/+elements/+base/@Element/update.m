function update(obj)
    %UPDATE (protected) Add tags to the elements and sides
    %   This tags element objects based on the tags contained in the nodes 
    %   associated with the element.
    %
    % Syntax
    %   update()
    %
    % Description
    %   update() updates the tags for the element and side based on the
    %   node tags.
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

    % Loop through the elements
    for i = 1:length(obj);
        
        % Update element tags with all unique tags from the nodes
        tag = [obj(i).nodes.tag];
        [~,ix] = unique({tag.name});
        obj(i).tag = tag(ix);

        % Update sides
        updateSides(obj(i));
%         % Apply tag to the element
%         obj(i).tag{end+1} = tag;
% 
%         % Tag the sides
%         s_id = 1:obj(i).n_sides;
%         if boundary_flag;
%             s_id = s_id([obj(i).sides.on_boundary]);
%         end
% 
%         % Loop through sides, mark side if dofs match
%         for s = 1:length(s_id);
%             obj(i).sides(s_id(s)).tag{end+1} = tag;
%         end    
    end
end

function updateSides(elem)

    % Update the sides
    for s = 1:elem.n_sides;
       no = elem.nodes(elem.side_ids(s,:));
       
       % 1D case, single nodes on a side
       if length(no) == 1;
           elem.sides(s).tag = no.tag;
       else
            T1 = no(1).tag;         % tags from node 1 of side s
            T2 = [no(2:end).tag];   % tags from other nodes of side s
            
            % If the tag sets are not empty continue
            if ~isempty(T1) && ~isempty(T2);
                
                % Compare the name and type
                [Lia1,~] = ismember({T1.name},{T2.name});
                [Lia2,Locb] = ismember({T1.type},{T2.type});
                
                % Extract index of matching tags
                Lia = all([Lia1,Lia2],2);
                index = Locb(Lia);

                % tag the side if things match
                if ~isempty(index)
                    T0 = T2(index);
                    for i = 1:length(T0); % loop over all tags
                        if (strcmpi(T0(i).type,'boundary') && elem.sides(s).on_boundary) ||...
                           (strcmpi(T0(i).type,'subdomain') && ~elem.sides(s).on_boundary);
                                elem.sides(s).tag(end+1) = T0(i);
                        end
                    end
                end
            end
       end
    end
end

