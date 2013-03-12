function varargout = getPosition(obj,xi,varargin)
    %GETPOSITION Returns the real coordinates given xi, eta, ...
    %
    % Syntax
    %   x = ggetPosition(xi)
    %   [x,y,z] = getPosition(xi)
    %   [...] = getPosition(..., 'PropertyName', PropertyValue)
    %   xyz = getPosition(...)
    %
    % Description
    %   [...] = getPosition(...) returns the position in global
    %   system given a value(s) for the local position, the number 
    %   of outputs varies according the number of spatial 
    %   dimensions.
    %
    %   xyz = getPosition(...) same as above but it returns the
    %   positions as a single array.
    %
    %   [...] = getPosition(..., 'PropertyName', PropertyValue)
    %   allow the limiting of which shape functions are used for
    %   the mapping, see EXAMPLE14b for an example.
    %
    % GETPOSITION Property Descriptions
    %   index
    %       logical array | numeric array
    %       An array that limits what shape functions are used for
    %       the mapping, this is useful elements that have mutiple
    %       dofs per node. The value of index should be defined
    %       such that the following is valid
    %           N = obj.shape(xi,...)
    %           N = N(opt.index)   
    %       Using this property overrides the default behavior,
    %       which is to use opt.index = 1:n_dof_node:end when a the
    %       no. of nodes is different from no. of shape functions,
    %       as is the case for the Beam element. For the Beam
    %       element there are 2 dofs per node and the 1 and 3 value
    %       for the shape functions may be used for maping the
    %       displacement, this behavior is the default. If you have
    %       an element that behaves otherwise, use the index
    %       property.
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
    
    % Gather the options
    opt.index = [];
    opt = gatherUserOptions(opt,varargin{:});

    % The number of spatial dimensions available, no. dofs per node
    n = obj.n_dim;
    n_dof_node = obj.nodes(1).n_dof;

    % Loop through the dimensions and return the desired position
    xyz = zeros(1,n);
    coord = obj.nodes.getCoord();
    for i = 1:n

       % Evaluate shape functions
       N = obj.shape(xi,'-scalar');

       % The index property was set, so limit the shape functions
       if ~isempty(opt.index);
           N = N(opt.index);

       % The no. of shape functions is not equal to no. of nodes,
       % so assume that the 1:n_dof_node:end shape functions are used to
       % map the position, this is the case for Beam elements.
       elseif length(N) > length(varargin);
           N  = N(1:n_dof_node:end);
       end

       % Map the position
       xyz(i) = N*coord(:,i);
    end

    % Reduce to single array if only a single output is given
    if nargout > 1;
        varargout = num2cell(xyz,1);
    else
        varargout{1} = xyz;
    end
end