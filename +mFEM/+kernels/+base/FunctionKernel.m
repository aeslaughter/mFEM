classdef FunctionKernel < mFEM.kernels.base.Kernel;
%     properties(Access = public);
%         func;
%     end

    methods 
        function obj = FunctionKernel(name, input)
            obj = obj@mFEM.kernels.base.Kernel(name);
            
            if ~isa(input,'function_handle') && ~ischar(input);
                error('FunctionKernel:FunctionKernel', 'The input must be a function handle or a character string convertable to a handle, but a %s was given.',class(input));    
            end
            
            obj.value = input;
            
            if ischar(input);
                obj.value = str2func(['@(elem,x,t) ',input]);
            else
                obj.value = input;
            end
             
            n = nargin(obj.value);
            if n ~= 3;
                error('FunctionKernel:FunctionKernel', 'The function must take three variables (elem, x, t), the supplied function accepts %d.', n);    
            end
             
        end        

        function value = eval(obj,elem,x,varargin)
            if isempty(varargin); 
                t = [];
            else
                t = varargin{1};
            end
            
%             if ischar(obj.func);
%                 fcn = str2func(['@(elem,x,t) ',obj.func]);
%             else
%                 fcn = obj.func;
%             end
%              
%             n = nargin(fcn);
%             if n ~= 3;
%                 error('FunctionKernel:FunctionKernel', 'The function must take three variables (elem, x, t), the supplied function accepts %d.', n);    
%             end
            
            value = feval(obj.value,elem,x,t);
        end

        
    end

end