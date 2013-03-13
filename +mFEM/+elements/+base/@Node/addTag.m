function addTag(obj,tag)
    %ADDTAG Appends the supplied tag to the node tags
    %   The Mesh class allows for boundaries and subdomains to be tagged,
    %   this is stored on each node to control integration for finite
    %   element solutions.
    %
    % Syntax
    %   addTag(tag)
    %   
    % Description
    %   addTag(tag) appends the supplied tag to a cell array of tags.
    %
    % See Also Mesh
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
    for i = 1:length(obj);
       obj(i).tag(end+1) = tag;
    end
end