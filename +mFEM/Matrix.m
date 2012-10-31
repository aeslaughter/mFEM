classdef Matrix < mFEM.handle_hide

   properties (Access = private)
      I;
      J;
      Aij;
      m
      n
   end
   
   methods
       function obj = Matrix(varargin)
           %MATRIX Class constructor
           %
           % Syntax:
           %    Matrix(mesh)
           %    Matrix(m,n)
           %    Matrix(m)
           
           if nargin == 1 && isa(varargin{1},'mFEM.FEmesh');
               obj.m = varargin{1}.n_dof;
               obj.n = obj.m;
           elseif nargin == 1;
               obj.m = varargin{1};
               obj.n = obj.m;
           elseif nargin == 2;
               obj.m = varargin{1};
               obj.n = varargin{2};
           else
               error('Matrix:Matrix', 'Input Error');
           end
       end
       
       function varargout = size(obj, varargin)
           %SIZE
           %
           % Syntax:
           %    [m,n] = size()
           %    m = size()
           %    m = size(dim)

           out = [obj.m, obj.n];
           
           if nargin == 1 && nargout <= 1;
               varargout{1} = out;

           elseif nargin == 1 && nargout == 2;
               varargout{1} = out(1);
               varargout{2} = out(2);
               
           elseif nargin == 2 && nargout <= 1;
               varargout{1} = out(varargin{1});
   
           else
               error('Matrix:size','Input error');
           end
       end
       
       function add_matrix(obj, B, varargin)
           %ADD_MATRIX Adds a sub-matrix to the global matrix
           %
           % Syntax
           %    add_matrix(B,dof)
           %    add_matrix(B,dof1,dof2)
           
           % Parse the input
           
            if nargin == 3;
                dof1 = varargin{1};
                dof2 = dof1;
            elseif nargin == 4;
                dof1 = varargin{1};
                dof2 = varargin{2};
            end

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
           A = sparse(obj.I, obj.J, obj.Aij);
       end

   end
end