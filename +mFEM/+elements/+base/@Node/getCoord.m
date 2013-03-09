function out = getCoord(obj)
    %GETCOORD Get the physical coordinates for the node
    %   The spatial coordinates for the node are stored as column vector in
    %   the node object, this function extracts the values. It functions on
    %   individual nodes or node arrays.
    %
    % Syntax
    %   out = getCoord()
    %
    % Description
    %   out = getCoord() retreives all of the spatial coordinates for each
    %   of the node objects and returns a numeric array, each columm
    %   representing the values for the corresponding node.
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
    out = [obj.coord];
end