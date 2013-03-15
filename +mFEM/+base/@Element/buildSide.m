function side = buildSide(obj,id)
    %BUILDSIDE Build an element for the side
    %
    % Syntax
    %   side = buildSide(id)
    % 
    % Description
    %   side = buildSide(id) creates an element for the specified
    %   side, the type of element is specified in the side_type
    %   property.
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
    
    % Extract the nodes for the side
    ix = obj.side_ids(id,:);
    node = obj.nodes(ix);

    % Create the side element
    side = feval(['mFEM.elements.',obj.side_type],NaN,node);

    % Set the global dofs
    side.dof = [side.nodes.dof];
end