classdef Matrix < handle

   properties %(Access = private)
      I;
      J;
      Aij;
      m;
      n;
      
   end
   
   methods
       function obj = Matrix(varargin)
           if nargin == 1;
               obj.m = varargin{1};
               obj.n = obj.m;
           elseif nargin == 2;
               obj.m = varargin{1};
               obj.n = varargin{2};
           else
               disp('error');
           end
%            

           
         
       end
       
%        function varargout = size(obj, dim)
% 
%            if nargin == 1 && nargout <= 1;
%                if nargout == 0;
%                    [obj.m, obj.n]
%                else
%                     varargout{1} = [obj.m, obj.n];
%                end
%                
%            elseif nargin == 1 && nargout == 2;
%                varargout{1} = obj.m;
%                varargout{2} = obj.n;
%                
%            elseif nargin == 2 && nargout <= 1;
%                 out = [obj.m, obj.n];
%                 
%                 if nargout == 0;
%                     out(dim)
%                 else
%                     varargout{1} = out(dim);
%                 end
%            else
%                disp('error...');
%            end
%        end
       
       function add_matrix(obj, B, varargin)
           
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
            obj.I(idx) = repmat(dof1, length(dof1),1);
            obj.J(idx) = sort(repmat(dof2, length(dof2), 1));

            % Add the local mass and stiffness matrix to the sparse matrix values
            obj.Aij(idx) = reshape(B, numel(B), 1);

       end
       
       function A = build(obj)

            	A = sparse(obj.I, obj.J, obj.Aij);

       end
       
       

       
   end
    
    
    
    
    
    
end