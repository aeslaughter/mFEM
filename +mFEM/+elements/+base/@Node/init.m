function init(obj,id,x,varargin) 
    %INIT Initializes the node objects, used for creating arrays quickly
    %   Allows for the objects in array of objects to be looped over
    %   internally, this is a significantly faster (order of magnitude)
    %   method for initializing classes.
    %
    % Syntax
    %   init(id,nodes)
    %   init(id,nodes,space)
    %
    % Description
    %   init(id,nodes) loops over the node objects and sets the ids to the
    %   values in the id array and coordinate array of nodes. The length of
    %   id and the number of rows in nodes must be equal to the number of
    %   node objects.
    %
    %   init(id,nodes,space) also allows for the number of
    %   degrees-of-freedom per node to be set. The value of space may be a
    %   scalar (i.e., 1,2,3,...) indicating the number of dofs for the node
    %   or it may be the keyword 'vector' indicating the the no. of dofs is
    %   equal to the no. of space dimensions, which is the no. of columns
    %   of the nodes input.
    %
    % See Also Element
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
    
    % Get the dimensions for the node array
    [n,ndim] = size(x);
    
    % Check that the correct no. of node objects exists
    if n ~= length(obj);
        error('Node:init:InvalidCoordDimension','The node inputs are not sized correctly');
    end
    
    % Check that the no. of ids is correct
    if n ~= length(id);
        error('Node:init:InvalidIDDimension','The id input is not sized correctly.');
    end
    
    % Handle the numeric space input
    if nargin == 4 && isnumeric(varargin{1});
       ndof = varargin{1};
       
    % Handle the keyword input
    elseif nargin == 4 && ischar(varargin{1}) && strcmpi('vector',varargin{1});
       ndof = ndim;
       
    % The default case
    else
       ndof = 1;
    end

    % Setting the properties of the nodes is faster when set using the
    % bracket notation avaiable for structs and classes. This requires that
    % the input be put into a cell.
%     x = mat2cell(x',ndim,ones(1,n));        % node coordinates
    n_dim = num2cell(repmat(ndim,n,1));      % no. of space dim.
    n_dof = num2cell(repmat(ndof,n,1));      % no. of dofs per node
    proc = num2cell(repmat(labindex,n,1));  % the processor
    id = num2cell(id);                      % global id

    % Pass the input to the objects
    [obj.id] = id{:};
    [obj.n_dim] = n_dim{:};
    [obj.n_dof] = n_dof{:};
    [obj.lab] = proc{:};
    %[obj.coord] = x{:};    
    
    % Loop through nodes and assign coordinates (req. for padding zeros)
    for i = 1:n;
        obj(i).coord(1:ndim) = x(i,:);
    end
end