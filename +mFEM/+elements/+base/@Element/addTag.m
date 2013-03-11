function addTag(obj,tag,varargin)
    %ADDTAG (protected) Add tags to the elements and sides
    %   This adds tags element objects, this function was designed to be
    %   used by the mFEM.Mesh classes own addTag method. Use it directly
    %   with caution, especially with the 'boundary' flag which assumes
    %   only elements on the boundary are supplied. Thus, the supplied tag
    %   is always attached to the element without checking if the element
    %   is on the boundary. This was done for speed reasons.
    %
    % Syntax
    %   addTag(tag)
    %   addTag(tag,'boundary')
    %
    % Description
    %   addTag(tag) adds the tag (char) to all of the elements and sides
    %   addTag(tag,'boundary') adds the tag to all of the elements but
    %   limits the application to sides on the boundary, this version is
    %   used by the 
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
        
    % Parse the input
    boundary_flag = false;
    if nargin == 3 && strcmpi(varargin{1},'boundary');
        boundary_flag = true;
    end
    
    % Loop through the elements
    for i = 1:length(obj);

        % Apply tag to the element
        obj(i).tag{end+1} = tag;

        % Tag the sides
        s_id = 1:obj(i).n_sides;
        if boundary_flag;
            s_id = s_id([obj(i).sides.on_boundary]);
        end

        % Loop through sides, mark side if dofs match
        for s = 1:length(s_id);
            obj(i).sides(s_id(s)).tag{end+1} = tag;
        end    
    end
end