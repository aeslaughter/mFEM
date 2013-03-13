function tagEmptyBoundary(obj, tag)
    %TAGEMPTYBOUNDARY Tag all nodes on boundary not previously tagged
    %
    % Syntax
    %   tagEmptyBoundary(id)
    %
    % Description
    %   tagEmptyBoundary(id) tags all the nodes on boundaries that are 
    %   unmarked with the specified id.
    %
    % See Also ADDBOUNDARY
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
    
    % Convert the tag to a string, if it is not
    if isnumeric(tag); 
        tag = num2str(tag);
    end
    
    % Build tag structure
    tag = struct('name',tag,'type','boundary');
    
    % spmd statements do not accept dot referencing
    nodes = obj.nodes;
    
    % Begin parralel computations
    spmd
        % Loop through all local elements
        for i = 1:length(nodes);
            
            % The current element
            no = nodes(i);
            
            % Go the the next node if it is not on a boundary or if it
            % already has a tag applied
            if ~no.on_boundary || ~isempty(no.tag);
                continue;
            end
            
            % Apply the tag to the element
            no.tag = tag;
        end
    end
    
    % Return the node objects
    obj.nodes = nodes;
end