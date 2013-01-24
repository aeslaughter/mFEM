classdef Constant < mFEM.kernels.base.Kernel;
 
    properties (Access = public)%(Access = ?mFEM.registry.base.Registry)
        value;
        reserved = {};
    end
    
    methods 
        function obj = Constant(name, input)
            obj = obj@mFEM.kernels.base.Kernel(name);
            obj.testName(name);
            if isnumeric(input);
                obj.value = mat2str(input);
            elseif ischar(input)
                obj.value = input;
            else
                error('Constant:Constant', 'The input must be a numeric or character, but a %s was given.', class(input));    
            end
        end        

        function str = apply(obj, str)
            expr = ['\<',obj.name,'\>'];
            repstr = obj.value;

            if ~ischar(str);
                error('Kernel:apply', 'The supplied input (str) must be a character string');
            end 

            str = regexprep(str, expr, repstr); 
        end
        
        function value = get(obj)
            value = obj.eval();
        end
        
        function value = eval(obj,varargin)
            value = eval(obj.value);
        end  
    end
end