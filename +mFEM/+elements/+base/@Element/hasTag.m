function varargout = hasTag(obj,tag)
    %HASTAG Returns true if element as the supplied tag
    %
    % Syntax
    %   tf = hasTag(tag)
    %   [tf,sid] = hasTag(...)
    %
    % Description
    %   tf = hasTag(tag) returns a true if the element contains the given
    %   tag, which was added with addBoundary or addSubdomain.
    %
    %   [tf,sid] = hasTag(...) same as above but also returns a vector of
    %   integers indicating which side(s) also have the tag.
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

    % Convert numerics to char, throw error if it is something else
    if isnumeric(tag); 
        tag = num2str(tag); 
    elseif ~ischar(tag);
        error('Element:hadTag:InavlidInput','The tag must be a number of character');
    end

    % Single output case, return true if the tag is found
    if nargout == 1;
        varargout{1} = any(strcmp(tag,obj.tag));

    % Dual output case
    else
        % Initialize tf and side index storage
        tf = false(obj.n_sides,1);
        sid = 1:obj.n_sides;
        
        % Loop through the sides, searching for tag
        for i = 1:obj.n_sides
            tf(i) = any(strcmp(tag,obj.sides(i).tag));
        end
        
        % Return the results
        varargout{1} = any(tf);
        varargout{2} = sid(tf);
    end
end