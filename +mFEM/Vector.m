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
    
    properties %(Access = private)
        f     % the vector
        codist;
        is_parallel;
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
               obj.m = m.n_dof;
               obj.initialize();
               obj.zero();
               
           % Case when only m is specfied     
           elseif nargin == 1 && isscalar(m);
               obj.m = m;
               obj.initialize();
               obj.zero();

           % Case when a vector is given    
%            elseif nargin == 1 && iscodistributed(m);
%                
%                obj.m = length(m);
%                obj.f = m;
               
               
               
           % Case when a constant vector is desired
           elseif nargin == 2;
               obj.m = m;
               obj.initialize();
               obj.setConstant(varargin{1});

           % Input not understood
           else
               error('Vector:Vector', 'Input Error');
           end
           


       end
       
       function zero(obj)
           
           if obj.is_parallel;
                codist = obj.codist;
                m = obj.m;
                spmd
                    f = codistributed.zeros([m,1], codist);
                end
                obj.f = f;
           else
              obj.f = zeros(obj.m,1); 
           end
 
       end
       
       function value = get(obj,in)
           %GET Extract values from the vector
           %
           % Syntax
           %    x = get(idx)
        
        n = length(in);   
           
        if any(in < 1) || any(in > obj.m);
            error('Vector:get:OutOfRange', 'The supplied index value(s) must in range of the data (1 to %d)',obj.m);
        end
           
        if obj.is_parallel;
            vec = obj.f;
            spmd  
                
                idx = globalIndices(vec,1);
                
                pos = NaN(n,1);
                for i = 1:length(in);
                    p = find(idx == in(i));
                    if ~isempty(p); pos(i) = p; end
                
                end
                gsize = [n*numlabs,1];


                codist = codistributor1d(codistributor1d.unsetDimension, ...
                    codistributor1d.unsetPartition, gsize);
                P = codistributed.build(pos, codist); 
            end

            value = gather(vec(~isnan(P)));
        else
            value = obj.f(in);
        end
        
       end
       
       function init(obj,local)
          %INIT
%           
%           f = obj.f
%           spmd
%             codist = getCodistributor(f) 
%             f = codistributed.build(local, codist);
%           end
% %           obj.f = f;
%           
       end
    end
    
    methods (Access = private)
        
        function initialize(obj)
           
           if matlabpool('size') > 1;
               obj.is_parallel = true;
            
               m = obj.m;
               spmd
                   codist = codistributor1d(codistributor1d.unsetDimension, ...
                                codistributor1d.unsetPartition, [m,1]);
               end
               obj.codist = codist;
               
           else
               obj.is_parallel = false; 
               obj.codist = [];
           end
        end
        
        function setConstant(obj, c)
            
            if obj.is_parallel;
                codist = obj.codist;
                m = obj.m;

                spmd
                    n = codist.Partition(labindex);
                    local = ones(n,1)*c;
                    f = codistributed.build(local, codist); 
                end
                obj.f = f;
            else
               obj.f = ones(obj.m,1)*c; 
            end
            
        end

   end
end