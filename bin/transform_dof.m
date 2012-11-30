function D = transform_dof(d, n)
    %TRANSFROM_DOF Converts the dofs for vector element space
    %
    % Syntax
    %   D = transform_dof(d,n)
    %   
    % Description
    %   D = transform_dof(d,n) converts the scalar degrees of freedom
    %       for to vector based degrees of freedom. For example,
    %       inputing d = [1,3], n = 2 returns D = [1,2,5,6].
    %
    %----------------------------------------------------------------------
    %  mFEM: An Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2012 Andrew E Slaughter
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
    %  along with this program.  If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
    
    % Size of vector space dofs
    D = zeros(n*length(d),1); 

    % Loop through dimensions and build vector
    for i = 1:n;
        D(i:n:end) = d*n - (n-i);
    end 