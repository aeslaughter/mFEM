function out = getParents(obj,varargin)
    %GETPARENTS Gets the element handles listed as parents to the node(s)
    %   The Mesh class uses this feature when creating elements, it allows
    %   for the neighbor finding method (see Element) to focus efforts
    %   making the process much faster.
    %
    % Syntax
    %   getParent()
    %   getParent(elem)
    %
    % Description
    %   getParent() get the complete list of parent elements
    %
    %   getParent(elem) gets the list of parent elements, excluding the
    %   element given.
    %
    % See Also Element addParents
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

    % Get the list of parents
    out = [obj.parents];
    
    % Return if it is empty
    if isempty(out); 
        return
    end
    
    % Get the unique values
    [~,ix] = unique([out.id]);
    out = out(ix);

    % Exclude the supplied element, if present
    if nargin == 2;
        out = out(out~=varargin{1});
    end
end