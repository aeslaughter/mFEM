function N = shape(obj, x, varargin)
    %SHAPE Returns the shape functions
    %
    % Syntax
    %   shape(xi)
    %   shape(xi, '-scalar')
    %
    % Description
    %   shape(xi) returns the element shape functions evaluated at
    %   the local locations specified by xi.
    %
    %   shape(xi,'-scalar') allows user to
    %   override the vectorized output using the scalar flag, this
    %   is used by GETPOSITION
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
    
    % Determine no. of dofs per node
    n_dof_node = obj.nodes(1).n_dof; % no. of dofs per node

    % Parse options (do not use gatherUserOptions for speed)
    scalar_flag = false;
    if nargin == 3 && strcmpi(varargin{1},'-scalar');
        scalar_flag = true;            
    end                

    % Scalar field basis functions
    N = obj.basis(x);

    % Non-scalar fields
    if ~scalar_flag && n_dof_node > 1;
        n = N;                          % re-assign scalar basis
        r = n_dof_node;                 % no. of rows
        c = n_dof_node*obj.n_nodes;     % no. of cols
        N = zeros(r,c);                 % size the vector basis

        % Loop through the rows and assign scalar basis
        for i = 1:r;
            N(i,i:r:c) = n;
        end
    end      
end