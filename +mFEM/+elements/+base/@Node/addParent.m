 function addParent(obj,elem)
    %ADDPARENT Adds supplied element handle to list of parent elements
    %   The Mesh class uses this feature when creating elements, it allows
    %   for the neighbor finding method (see Element) to focus efforts
    %   making the process much faster.
    %
    % Syntax
    %   addParent(elem)
    %
    % Description
    %   addParent(elem) adds the supplied element to a list of elements
    %
    % See Also Element getParents
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
    
    % Loop through the nodes and create or append the element to the list
    for i = 1:length(obj);
        idx = length(obj(i).parents);
        if idx == 0;
            obj(i).parents = elem;
        else
            obj(i).parents(idx+1) = elem;
        end
    end
end 
        
