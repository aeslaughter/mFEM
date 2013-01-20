classdef Matrix < handle
    %MATRIX A wrapper class for MATLAB's sparse matrix creation
    % It is best to create the sparse matrices in MATLAB using an index 
    % based assembly, see doc sparse. This class automates the generation
    % of the i,j,s vector creation. The main purpose is to provide a means
    % of adding a dense matrix into the sparse based on the degrees of
    % freedom, i.e. inserting the element matrix into the global.
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
      I;    % i vector (see doc sparse)
      J;    % j vector (see doc sparse)
      Aij;  % s vector (see doc sparse)
      m     % no. for rows
      n;    % no. of columns
    end
   
    methods
       function obj = Matrix(varargin)
           %MATRIX Class constructor
           %
           % Syntax
           %    Matrix(mesh)
           %    Matrix(m,n)
           %    Matrix(m)
           %
           % Description
           %    Matrix(mesh) create a global matrix based on the mesh
           %    contained in the FEMESH object.
           %
           %    Matrix(m,n) create an m x n matrix.
           %    
           %    Matrix(m) create an m x m matrix.

           % Case when creating from an FEmesh object
           if nargin == 1 && isa(varargin{1},'mFEM.FEmesh');
               mesh = varargin{1};
               obj.m = mesh.n_dof;
               obj.n = obj.m;
               
%                N = 0;
%                for e = 1:mesh.n_elements;
%                    N = N + mesh.element(e).n_dof^2;
%                end
%                obj.I = zeros(N,1);
%                obj.J = zeros(N,1);
%                obj.Aij = zeros(N,1);
           
           % Case when only m is specfied     
           elseif nargin == 1;
               obj.m = varargin{1};
               obj.n = obj.m;
               
           % Case when both m and n are given    
           elseif nargin == 2;
               obj.m = varargin{1};
               obj.n = varargin{2};
           else
               error('Matrix:Matrix', 'Input Error');
           end
       end
       
       function varargout = size(obj, varargin)
           %SIZE Returns the size of the the matrix
           %
           % Syntax
           %    [m,n] = size()
           %    m = size()
           %    m = size(dim)
           %
           % Description
           %    [m,n] = size() returns the size of the matrix as an
           %    independant outputs.
           %
           %    m = size() returns the size as a vector, m = [m,n]
           %
           %    m = size(dim) returns the dim compontent of the dimensions,
           %    dim = 1 returns m, dim = 2 returns n.

           % The size of the matrix
           out = [obj.m, obj.n];
           
           % Output in vector format
           if nargin == 1 && nargout <= 1;
               varargout{1} = out;

           % Output two variables (m,n)
           elseif nargin == 1 && nargout == 2;
               varargout{1} = out(1);
               varargout{2} = out(2);
               
           % Output the dimension of interest
           elseif nargin == 2 && nargout <= 1;
               varargout{1} = out(varargin{1});
   
           else
               error('Matrix:size','Input error.');
           end
       end
       
       function add(obj, B, varargin)
            %ADD Adds a dense sub-matrix to the global matrix
            %
            % Syntax
            %    add(B, dof)
            %    add(B, dof1, dof2)
            %
            % Description
            %    add(B, dof) adds the dense matrix B to the locations
            %    specified in dof, this is equivelent to the following:
            %        M(dof,dof) = M(dof,dof) + B,
            %    where M is the global sparse matrix.
            %
            %    add(B, dof1, dof2) adds the dense matrix B to the 
            %    locations specified in dof1 and dof2, this is equivelent to 
            %    the following:
            %        M(dof1,dof2) = M(dof1,dof2) + B,
            %    where M is the global sparse matrix.
           
            % Case when only a single dof vector is supplied
            if nargin == 3;
                dof1 = varargin{1};
                dof2 = dof1;
                
            % Case when both dof vectors are given    
            elseif nargin == 4;
                dof1 = varargin{1};
                dof2 = varargin{2};
            end

            % Make sure the dof vectors are organized as columns
            if ~iscolumn(dof1); dof1 = dof1'; end
            if ~iscolumn(dof2); dof2 = dof2'; end

            % Compute indices for inserting into sparse matrix i,j,s vectors
            l = length(obj.I);
            b = numel(B);
            idx = l+1 : l+b;

            % Build the i,j components for the sparse matrix creation
            i = repmat((1:length(dof1))', length(dof1), 1);
            j = sort(i);
            obj.I(idx) = dof1(i);
            obj.J(idx) = dof2(j);

            % Add the local mass and stiffness matrix to the sparse matrix values
            obj.Aij(idx) = reshape(B, numel(B), 1);
       end
       
       function A = init(obj)
           %INIT Create and return the sparse matrix
           %
           % Syntax
           %    init()
           %
           % Description
           %    init() creates and returns the sparse matrix from the 
           %    I,J,and Aij vectors.
           A = sparse(obj.I, obj.J, obj.Aij);
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
          
           % Delete I,J,Aij properties
           obj.I = []; obj.J = []; obj.Aij = [];
       end
   end
end