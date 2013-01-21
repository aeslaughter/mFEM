classdef FunctionKernel < mFEM.kernels.base.Kernel;
    properties
        input;  % function as input
        constReg;
        options = struct('constants',[]);
    end

    methods 
        function obj = FunctionKernel(name, input, varargin)
            obj = obj@mFEM.kernels.base.Kernel(name);
              
            obj.options = gather_user_options(obj.options,varargin{:});
            if isa(obj.options.constants, 'mFEM.registry.ConstantKernelRegistry');
               obj.constReg = obj.options.constants;
            end
            
            if ~isa(input,'function_handle') && ~ischar(input);
                error('FunctionKernel:FunctionKernel', 'The input must be a function handle or a character string convertable to a handle, but a %s was given.', class(input));    
            end
            
            obj.input = input;
            
%             if ischar(obj.func);
%                 if ~isempty(obj.constReg)
%                     obj.constReg.apply(obj);
%                 end
%                 obj.func = str2func(['@(elem,x,t) ',obj.func]);
%             else
%                 obj.func = input;
%             end
%              
%             n = nargin(obj.func);
%             if n ~= 3;
%                 error('FunctionKernel:FunctionKernel', 'The function must take three variables (elem, x, t), the supplied function accepts %d.', n);    
%             end
        end        
        
        function str = apply(obj, str, elem, x, t)
            
            if ~ischar(str);
                error('FunctionKernel:apply', 'The input must be a character string');
            end    
            
            % Apply OBJ's value to KERN
            expr = ['\<',obj.name,'\>'];
            repstr = obj.eval(elem, x, t);
            
            str = regexprep(str, expr, repstr); 
        end
        
        function value = eval(obj,elem,x,varargin)
            if isempty(varargin); 
                t = [];
            else
                t = varargin{1};
            end
            
            obj.value = obj.input;
            if ischar(obj.value);
                if ~isempty(obj.constReg)
                    obj.constReg.apply(obj);
                end
                fcn = str2func(['@(elem,x,t) ',obj.value]);
            else
                fcn = obj.value;
            end
             
            n = nargin(fcn);
            if n ~= 3;
                error('FunctionKernel:FunctionKernel', 'The function must take three variables (elem, x, t), the supplied function accepts %d.', n);    
            end
   
            value = feval(fcn, elem, x, t);
        end

        
    end

end