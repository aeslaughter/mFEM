function setBoundaryFlag(obj)
    %SETBOUNDARYFLAG Sets the boundary flag to true for the elements
    %   Each node has a true|false flag indicating if it is on a boundary,
    %   which is defined whether or not it a part of a shared side. This
    %   flag is set by the Element class.
    %
    % Syntax
    %   setBoundarFlag()
    %
    % Description
    %   setBoundaryFlag() sets the on_boundary flag for the elements
    %
    % See Also Element Mesh
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
    t = num2cell(true(length(obj),1));    
    [obj.on_boundary] = t;
end