function init(obj,id,nodes,varargin)
    %INIT Element initialization for creating vectors of element objects
    %
    % Syntax
    %   init(id,nodes)
    %   init(id,nodes,'-addParents')
    %
    % Description
    %   init(id,nodes) creates an element object for each id and set of
    %   Node objects. Each element is constructed from the rows of nodes,
    %   which should be the length of the ids.
    %
    %   init(id,nodes,'-addParents') same as above but triggers a call to
    %   the addParents method of the nodes, this is used by mFEM.Mesh to
    %   speed up neighbor finding.
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
    
    % The no. of elements
    n = length(obj);

    % Set/find the known constants and quadrature information
    lab = labindex; 
    
    if ~isempty(obj(1).quad);
        [qp,W] = obj(1).quad.rules('-cell');
        qp = repmat({qp},n,1);
        W = repmat({W},n,1);
        [obj.qp] = qp{:};
        [obj.W] = W{:};
    end
    
    % Accunt for single node input, make sure it is column
    if size(nodes,1) == obj(1).n_nodes && size(nodes,2) == 1;
        nodes = nodes';
    end
    
    % Build correctly sized cell array
    n_sides = size(obj(1).side_ids,1);
    n_sides = num2cell(repmat(n_sides,n,1));
    lab = repmat({lab},n,1);
    id = num2cell(id);

    % Set properties of elements
    [obj.n_sides] = n_sides{:};
    [obj.id] = id{:};
    [obj.lab] = lab{:};  

    % Loop through elements and set remaining properties (todo)
    for i = 1:length(obj);
        obj(i).nodes = nodes(i,:);
        obj(i).sides = struct('neighbor',[],...
                         'neighbor_side',[],...
                         'on_boundary',...
                         num2cell(true(obj(i).n_sides,1)),...
                         'tag',[]);
        obj(i).dof = [obj(i).nodes.dof];
        obj(i).n_dof = sum([obj(i).nodes.n_dof]);
    end
    
    if ~isempty(varargin) && strcmpi(varargin{1},'-addParents');
        for i = 1:length(obj);
            obj(i).nodes.addParent(obj(i));
        end   
    end
end