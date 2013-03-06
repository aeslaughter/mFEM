function tagEmptyBoundary(obj, tag)
    %TAGEMPTYBOUNDARY Mark all unmarked boundaries
    %
    % Syntax
    %   tagEmptyBoundary(id)
    %
    % Description
    %   tagEmptyBoundary(id) tags all boundaries that are unmarked
    %   with the specified id.
    %
    % See Also ADDBOUNDARY
    %
    %----------------------------------------------------------------------
    %  mFEM: An Object-Oriented MATLAB Finite Element Library
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
    %  along with this program.  If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
    
    % Convert the tag to a string, if it is not
    if isnumeric(tag); 
        tag = num2str(tag);
    end
    
    % spmd statements do not accept dot referencing
    elements = obj.elements;
    
    % Begin parralel computations
    spmd
        % Loop through all local elements
        for e = 1:length(elements);
            
            % The current element
            elem = elements(e);
            
            % Go the the next element if it is not on a boundary or if it
            % already has a tag applied
            if ~elem.on_boundary || ~isempty(elem.tag);
                continue;
            end
            
            % Apply the tag to the element
            elem.tag = tag;
            
            % Loop through the sides and apply tag
            for s = 1:elem.n_sides;
                
                % Go to the next side if not on the boundary
                if ~elem.sides(s).on_boundary;
                    continue;
                end
            
                % Apply the tag to the side
                elem.side(s).tag = tag;
            end
        end
    end
end