function B = shapeDeriv(obj,x)
    %SHAPEDERIV Returns shape function derivatives in global x,y system
    %
    % Syntax
    %   shapeDeriv(x)
    %
    % Description
    %   shapeDeriv(x) returns the element shape function 
    %   derivatives evaluated at the locations specified in xi.
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

    % Scalar field basis functin derivatives
    B = obj.gradBasis(x);

    % Non-scalar fields
    if n_dof_node > 1;
        b = B;                      % re-assign scalar basis derivatives
        r = obj.n_dof_node;         % no. of rows
        c = r*size(b,2);            % no. of cols
        B = zeros(r+1,c);           % size the vector basis

        % Loop through the rows and assign scalar basis
        for i = 1:r;
            B(i,i:r:c)  = b(i,:);
            B(r+1, i:r:c) = b((r+1)-i,:);
        end
    end
end