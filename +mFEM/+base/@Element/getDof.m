function dof = getDof(obj, varargin)
    %GETDOF The global degrees of freedom, account for type of space
    %
    % Syntax
    %   dof = getDof()
    %   dof = getDof('PropertyName', PropertyValue,...)
    % 
    % Description           
    %   dof = getDof() returns all of the global dofs for element
    %
    %   dof = getDof('PropertyName', PropertyValue) returns the
    %   dofs for the element subject to the restrictions specified
    %   by the properties, see the descriptions below for details.  
    %
    % GETDOF Property Descriptions
    %   Side
    %       integer
    %       Indicates that only the degrees of freedom for the
    %       specific side should be returned.
    %
    %   Local
    %       true | {false}
    %       Toggles the type of degrees of freedom to return, which
    %       is generally import for sides. For example, the
    %       following returns the local degrees of freedom for side
    %       number 1 of an element.
    %           dof = getDof('Side',1,'local',true) or
    %           dof = getDof('Side',1,'-local')
    %
    %   Component
    %       scalar | 'x' | 'y' | 'z'
    %       Returns the dof associated with the vector space or
    %       multipe dof per node elements like the Beam element.
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
    
    % Set default options and gather the options
    opt.side = [];
    opt.local = false;
    opt.component = [];
    opt = gatherUserOptions(opt, varargin{:});
    
    % Local only operates on single element
    if (~isempty(opt.side) || opt.local) && length(obj) > 1;
        error('Element:getDof:InvalidUseOfLocal','The local and side properties only functions on a single element');
    end
    
    % Determine the nodes to operate on
    if ~isempty(opt.side);
        nodes = obj.nodes(obj.side_ids(opt.side,:));
    else
        nodes = [obj.nodes];
    end
    
    % Extract the dofs
    dof = nodes.getDof('Component',opt.component); % desired dofs
    
    % Covert the global, to local if desried
    if opt.local
        glbl = obj.dof;            % all dofs, global
        ix = ismember(glbl,dof);   % location in glbl of desired dofs
        lcl = 1:length(obj.dof);   % local dofs
        dof = lcl(ix);             % desired local dofs
    end
end