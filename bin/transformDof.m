function D = transformDof(d, n)
    %TRANSFROMDOF Converts the dofs for vector element space
    %  Performs a simple mapping of a single degree-of-freedom (i.e.,
    %  the node id) to a multiple degree-of-freedom system.
    %
    % Syntax
    %   D = transformDof(d,n)
    %   
    % Description
    %   D = transformDof(d,n) converts the scalar degrees of freedom
    %       for to vector based degrees of freedom. For example,
    %       inputing d = [1,3], n = 2 returns D = [1,2,5,6].
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
    
    % Do nothing if n == 1
    if n == 1; 
        D = d; 
        return 
    end
    
    % Size of vector space dofs
    D = zeros(n*length(d),1); 

    % Loop through dimensions and build vector
    for i = 1:n;
        D(i:n:end) = d*n - (n-i);
    end 