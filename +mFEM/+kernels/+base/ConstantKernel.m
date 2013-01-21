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

        
        function apply(obj,kern)
            % Apply OBJ's value to KERN
            expr = ['\<',obj.name,'\>'];
            repstr = obj.value;
            str = kern.value;
            
            if ~ischar(str);
                error('ConstantKernel:apply', 'The supplied function must be a character string');
            end 
            
            kern.value = regexprep(str, expr, repstr); 
            
%             switch class(kern); 
%                 case 'mFEM.kernels.base.ConstantKernel';
%                     kern.value = regexprep(str, expr, repstr); 
%                     
%                 case 'mFEM.kernels.base.FunctionKernel';
%                     if ~ischar(str);
%                         error('ConstantKernel:apply', 'The supplied function must be a character string');
%                     end
%                     kern.value = regexprep(str, expr, repstr);      
%                     
%                 case 'mFEM.kernels.AutoKernel';
%                     kern.value = re
%                     
%                 otherwise
%                     error('ConstantKernel:apply', 'Application of constants to the %s classes is not yet supported.', class(kern));
%             end
        end
        
        function value = eval(obj,varargin)
            value = eval(obj.value);
        end

        
    end

end