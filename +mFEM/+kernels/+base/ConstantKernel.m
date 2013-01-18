classdef ConstantKernel < mFEM.kernels.base.Kernel;
    properties(Access = public);
       value;

    end

        
    methods 
        function obj = ConstantKernel(name, input)
            obj = obj@mFEM.kernels.base.Kernel(name,input);
            if isnumeric(input);
                obj.value = num2str(input);
            elseif ischar(input)
                obj.value = input;
            else
                error('ConstantKernel:ConstantKernel', 'The input must be a numeric or character, but a %s was given.', class(input));    
            end
        end        

        function merge(obj,kernel)
           
            if ~isa(kernel,class(obj));
               error('ConstantKernel:merge', 'Cannot merge with a class type of %s.', class(kernel));
            end
            
            obj.input = [num2str(obj.input),' + ', num2str(kernel.input)];
            obj.value = [obj.value, ' + ', kernel.value];
            obj.value = num2str(eval(obj.value));
            
        end
        
        function str = apply(obj,str)
            str = regexprep(str, ['\<',obj.name,'\>'], obj.value);  
        end
        
        function value = eval(obj,varargin)
            value = eval(obj.value);
        end

        
    end

end