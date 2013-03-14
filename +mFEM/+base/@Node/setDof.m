function dof = setDof(obj,varargin)
    %SETDOF (protected) sets the degrees-of-freedom for the node(s)
    %
    % Syntax
    %   setDof()
    %   setDof(strt)
    %   setDof(dof)
    %
    % Description
    %   setDof() sets one dof per node to the the global ids of the nodes
    %
    %   setDof(strt) sets the dofs for the nodes, in order of the nodes in
    %   the array starting with the value given.
    %
    %   setDof(dof) sets the dofs to the values to the values in the array,
    %   one row for each node.
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

    % Get the optional dof input
    if nargin == 2;
        dof = varargin{1};
    else
        dof = 1;
    end

    % Implicit case, increment from the starting value being sure to
    % account for nodes with varying dofs
    if isscalar(dof);
        strt = dof;
        for i = 1:length(obj);
            stop = strt + obj(i).n_dof - 1;
            obj(i).dof = uint32(strt:stop);
            strt = stop+1;
        end
    
    % Explicit case    
    elseif size(dof,1) == length(obj);
        for i = 1:length(obj);
            obj(i).dof = uint32(dof(i,:));
        end
    end
end