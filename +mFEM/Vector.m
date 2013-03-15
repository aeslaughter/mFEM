classdef Vector < mFEM.Matrix
    %VECTOR A class for creating a column vector.
    %   This limits the Matrix class to a column vector, it stores the
    %   vector as a sparse I,J,S format. This is necessary for this class
    %   to work properly with pVector (parallel Vector).
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
   
    methods
       function obj = Vector(m,varargin)
           %Vector Class constructor
           %
           % Syntax
           %    Vector(mesh)
           %    Vector(m)
           %
           % Description
           %    Vector(mesh) create a global vector based on the mesh
           %    contained in the FEMESH object.
           %
           %    Vector(m) create an m x 1 matrix.
           %
           %    Vector(m, c) creates an m x 1 matrix with all values
           %    assigned to c
           %
           %    Vector(vec) create an initilized vector

           % By default the I,J,S values of the sparse are empty
           I = [];
           J = [];
           Aij = [];
           
           % Case when creating from an Mesh object
           if nargin == 1 && isa(m,'mFEM.Mesh');
                m = m.n_dof;

           % Case when a vector is given    
           elseif nargin == 1 && ~isscalar(m);
                 [I,J,Aij] = find(m);
                 m = length(m);
                 
           % Case when a constant vector is desired
           elseif nargin == 2;
               [I,J,Aij] = find(ones(m,1)*varargin{1});
           end
           
           % Call the base constructor and set I,J,Aij values
           obj = obj@mFEM.Matrix(m,1);
           obj.I = I;
           obj.J = J;
           obj.Aij = Aij;
       end
       
       function add(obj,B,dof)
            %ADD Adds a sub-vector to the global vector
            %   This is an overloaded method of the matrix class, it simply
            %   limits the second dof input to 1 because this is a vector.
            %
            % Syntax
            %    add(fe,dof)
            %
            % Description
            %    add(fe,dof) adds the local vector fe to the 
            %    locations specified in dof, this is equivelent to the 
            %    following:
            %        f(dof) = f(dof) + fe,
            %    where f is the global vector.
            add@mFEM.Matrix(obj,B,dof,1);
       end
 
%        function out = getLocal(obj, dof)
%            %GETLOCAL extract a local vector given the dof
%            %
%            % Syntax
%            %    out = getLocal(dof)
%            %
%            % Description
%            %    out = getLocal(dof) returns a subvector, where dof are the
%            %    indices of the subvector
%            out = obj.f(dof);
%        end
   end
end