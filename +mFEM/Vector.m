classdef Vector < handle
    %VECTOR A wrapper class for creating a column vector, this is simple 
    % wrapper that provides the same features as the Matrix class so that 
    % syntax may remain consistent. This creates a column vector. And will
    % serve as the basis of creating parallel vectors in the future.  
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

    properties (SetAccess = private, GetAccess = public)
      m     % no. of rows
    end
    
    properties (Access = private)
        f     % the vector
    end
   
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

           % Case when creating from an FEmesh object
           if nargin == 1 && isa(m,'mFEM.FEmesh');
               mesh = m;
               obj.m = mesh.n_dof;
               obj.zero();
               
           % Case when only m is specfied     
           elseif nargin == 1 && isscalar(m);
               obj.m = m;
               obj.zero();
               
           % Case when a vector is given    
           elseif nargin == 1;
               obj.m = length(m);
               obj.f = m;
               
           % Case when a constant vector is desired
           elseif nargin == 2;
               obj.m = m;
               obj.f = ones(m,1)*varargin{1};

           % Input not understood
           else
               error('Vector:Vector', 'Input Error');
           end
       end
       
       function varargout = size(obj, varargin)
           %SIZE Returns the size of the the matrix
           %
           % Syntax
           %    m = size()
           %
           % Description
           %    m = size() returns the length of the vector

           % The size of the matrix
           varargout{1} = obj.m;
       end
       
       function add(obj, fe, dof)
            %ADD Adds a sub-vector to the global vector
            %
            % Syntax
            %    add(fe, dof)
            %
            % Description
            %    add(fe, dof) adds the local vector fe to the 
            %    locations specified in dof, this is equivelent to the 
            %    following:
            %        f(dof) = f(dof) + fe,
            %    where f is the global vector.

            obj.f(dof) = obj.f(dof) + fe;
       end
       
       function out = get_local(obj, dof)
           %GET_LOCAL extract a local vector given the dof
           %
           % Syntax
           %    out = get_local(dof)
           %
           % Description
           %    out = get_local(dof) returns a subvector, where dof are the
           %    indices of the subvector
           
           out = obj.f(dof);
       end
       
       function f = init(obj)
           %INIT Return the vector
           %
           % Syntax
           %    f = init()
           %
           % Description
           %    f = init() creturns the vector
           f = obj.f;
       end
       
       function zero(obj)
           %ZERO Clears the matrix
           %
           % Syntax
           %    zero()
           %
           % Description 
           %    zero() removes all existing values, it does not resize the
           %    matrix.

           obj.f = zeros(obj.m,1);
       end
   end
end