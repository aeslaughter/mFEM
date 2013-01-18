classdef ConstantKernel < mFEM.base.Kernel;
    properties(Access = public);
       input;
       value;

    end

        
    methods 
        function obj = ConstantKernel(name, input)
            obj = obj@mFEM.base.Kernel(name);
            if isnumeric(input);
                obj.input = num2str(input);
            elseif ischar(input)
                obj.input = input;
            else
                error('ConstantKernel:ConstantKernel', 'The input must be a numeric or character, but a %s was given.', class(input));    
            end
            obj.value = obj.input;
        end        

        function merge(obj,kernel)
           
            if ~isa(kernel,class(obj));
               error('ConstantKernel:merge', 'Cannot merge with a class type of %s.', class(kernel));
            end
            
            obj.input = [obj.input, ' + ', kernel.input];
            %obj.value = num2str(eval(obj.input));
            
        end
        
        function str = apply(obj,str)
            str = regexprep(str, ['\<',obj.name,'\>'], obj.input);  
        end
        
        function value = eval(obj,varargin)
            value = eval(obj.value);
        end

        
    end

end