classdef Constant < mFEM.kernels.base.Kernel;
 
    methods 
        function obj = Constant(name, input)
            obj = obj@mFEM.kernels.base.Kernel(name);
            if isnumeric(input) || ischar(input);
                obj.value = num2str(input);
            else
                error('Constant:Constant', 'The input must be a numeric or character, but a %s was given.', class(input));    
            end
        end        

        function value = eval(obj,varargin)
            value = eval(obj.value);
        end  
    end
end