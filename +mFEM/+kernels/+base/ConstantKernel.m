classdef ConstantKernel < mFEM.kernels.base.Kernel;
 
    methods 
        function obj = ConstantKernel(name, input)
            obj = obj@mFEM.kernels.base.Kernel(name);
            if isnumeric(input) || ischar(input);
                obj.value = num2str(input);
            else
                error('ConstantKernel:ConstantKernel', 'The input must be a numeric or character, but a %s was given.', class(input));    
            end
        end        

        function value = eval(obj,varargin)
            value = eval(obj.value);
        end  
    end
end