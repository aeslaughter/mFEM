function dof = getDof(obj,varargin)
    %GETDOF extract the global degrees of freedom for the node
    %   Each node holds an array of global degrees of freedom for the node,
    %   this function handles the node arrays when the nodes have varying
    %   spatial dimensions and degrees-of-freedom per node.
    %
    % Syntax
    %   dof = getDof()
    %   dof = getDof('PropertyName',PropertyValue,...)
    %
    % Description
    %   dof = getDof() returns a column vector containing all the dofs for
    %   tne nodes.
    %
    %   dof = getDof('PropertyName',PropertyValue,...) limits the dofs
    %   returned according to the criteria given.
    %
    % GETDOF Property Descriptions
    %   Component
    %       vector | scalar | 'x' | 'y' | 'z' | cell string
    %       Returns the dof associated with the vector space or
    %       multipe dof per node elements like the Beam element. Note, if
    %       you mesh was constructed with nodes that have varying dofs per
    %       node this property will cause errors.
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
    
    % Set the default options and gather those from the user
    opt.tag = {};
    opt.index = false;    
    opt.component = [];
    opt = gatherUserOptions(opt,varargin{:});
    
    % Get the component values
    cmp = obj.parseComponentInput(opt.component);
    
    % Case when all dofs are returned
    if isempty(cmp);
        dof = [obj.dof];  
  
    % Case when limited dofs are given
    else
         dof = {obj.dof};
         for i = 1:length(dof);
             dof{i} = dof{i}(cmp);
         end
         dof = cell2mat(dof); 
    end
end