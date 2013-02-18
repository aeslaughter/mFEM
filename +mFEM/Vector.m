classdef Vector < handle
    %VECTOR Creates a serial or parallel column vector
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

    properties (SetAccess = protected, GetAccess = public)
      m                         % no. of rows
      is_parallel = true;       % true if vector is codistributed
    end
    
    properties (Access = protected)
        f               % the vector (serial or codistributed)
        codist;         % the codistributor class used, if parallel 
    end
   
    methods
        function obj = Vector(m,varargin)
           %VECTOR Class constructor
           %
           % Syntax
           %    Vector(mesh)
           %    Vector(m)
           %    Vector(...,c)
           %    Vector(vec)
           %
           % Description
           %    Vector(mesh) create a global vector based on the mesh
           %    contained in the FEMESH object.
           %
           %    Vector(m) create an m x 1 matrix.
           %
           %    Vector(..., c) creates an m x 1 matrix with all values
           %    assigned to c
           %
           %    Vector(vec) create a vector from the numeric arrary vec

           if matlabpool('size') == 0;
              obj.is_parallel = false;
           end

           % Case when first argument is FEmesh object
           if nargin == 1 && isa(m, 'mFEM.FEmesh');
               obj.m = m.n_dof;     % Define the size from Mesh object
               obj.init(0);          % Initilize the vector

           % Case when first argument is a scalar (i.e., size)  
           elseif nargin >= 2 && isscalar(m) && isnumeric(varargin{1});
               obj.m = m;                % Define the size
               obj.init(varargin{1});    % Initilize

           % Case when a codistributed vector is given 
           elseif nargin == 1 && isdistributed(m); 
               vec = m;                 % Codistributed vector was given
               m = length(m);           % Define the length for use in spmd
               obj.m = m;               % Store the length in class

               % Begin parallel operations
               spmd
                   if ~iscolumn(vec); vec = vec'; end
                   codist = getCodistributor(vec);
               end

               % Store the codistributor and the distributed vector
               obj.codist = codist;
               obj.f = vec;
               
           % Case when serial vector is given   
           elseif nargin == 1 && ~iscodistributed(m);
               obj.m = length(m);   % define the vector length
               obj.init(m);         % initilize the vector             
               
           % Unknown input, produce an error
           else
               error('Vector:Vector:InvalidInput', 'The input for the Vector class constructor was not understood.');
           end
        end

        function zero(obj)
           %ZERO Sets all values of the vector to zero  
           %
           % Syntax
           %    zero();
           %
           % Description
           %    zero() sets all values of the vector to zero, it maintains
           %    the size and parallel distribution.
           obj.f(1:end) = 0; 
        end

        function value = get(obj,varargin)
            %GET Extract values from the vector
            %
            % Syntax
            %    x = get()
            %    x = get(idx)
            %
            % Description
            %    x = get() returns the vector, if the vector is serial this
            %    is simpy a numeric array. If the vector is parallel it
            %    returns the codistributed array.
            %
            %    x = get(idx) returns a numeric vector for the indices
            %    supplied in idx.
            
            % Return full vector
            if nargin == 1;
               value = obj.f;
               
            % Return specific values   
            else
               if any(varargin{1} < 1) || any(varargin{1} > obj.m);
                   error('Vector:get:OutOfRange', 'The index value(s) supplied are out of range for the vector.');
               end
               value = gather(obj.f(varargin{1}));
               if iscell(value) && isscalar(value);
                   value = value{1};
               end
            end
        end

        function add(obj, value, index)
           %ADD Adds subvector to the complete vector
           %
           % Syntax
           %    add(value, index)
           %
           % Description
           %    add(value, index) adds the value(s) to the global vector at
           %    the indices given in index. The value input may be a scalar
           %    or a numeric array the same size as index.
           obj.f(index) = obj.f(index) + value;
        end
    end
    
    methods (Access = private)     
        function init(obj, value)
            %INIT sets up the parallel/serial vector, initializing to zero
            %
            % Syntax
            %   init()
            %
            % Description
            %   init() creates a codistributor1d instance when
            %   running in parralel.
                                     
            % Build the vector
           if ~isscalar(value)
               if ~iscolumn(value); value = value'; end
           else
               value = repmat(value,[m,1]);
           end
                               
            % Parallel case
           % if matlabpool('size') > 0;
                obj.is_parallel = true; % set parallel flag to true

                % Begin parallel operations (create codistributor)
                spmd
                   vec = codistributed(value);
                   codist = getCodistributor(vec);
                end
                
                % Store the codistributor
                obj.codist = codist;
                obj.f = vec;
                
           % Serial case
%             else
%                obj.f = value;
%                obj.is_parallel = false; 
%                obj.codist = [];
%           end
        end
   end
end