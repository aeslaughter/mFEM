classdef ConstantKernel < mFEM.kernels.base.Kernel;
    
%     properties (Access = private)
%         value;
%     end
    
    methods 
        function obj = ConstantKernel(name, input)
            obj = obj@mFEM.kernels.base.Kernel(name);
            if isnumeric(input) || ischar(input);
                obj.value = num2str(input);
            else
                error('ConstantKernel:ConstantKernel', 'The input must be a numeric or character, but a %s was given.', class(input));    
            end
        end        

%         function merge(obj,kernel)
%            
%             if ~isa(kernel,class(obj));
%                error('ConstantKernel:merge', 'Cannot merge with a class type of %s.', class(kernel));
%             end
%             
%             obj.value = [obj.value,' + ', kernel.value];
%             obj.valur = num2str(eval(obj.value));
%             
%         end
        
%         function str = apply(obj,str)
%             % Apply OBJ's value to STR
%             if ischar(obj.value);
%                 str = regexprep(obj.value, ['\<',obj.name,'\>'], str); 
%             end
%         end
        
        function value = eval(obj,varargin)
            value = eval(obj.value);
        end

        
    end

end