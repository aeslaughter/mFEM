classdef ConstantKernel < mFEM.base.Kernel;
%     properties(Access = public);
%         value;
% 
%     end

        
    methods 
        function obj = ConstantKernel(system, name, input)
            obj = obj@mFEM.base.Kernel(system, name, input);
        end
        
        function output = init(obj)
            
            if isnumeric(obj.input);
                output = num2str(obj.input);
                return;
            end
            
            output = init@mFEM.base.Kernel(obj);
            output = num2str(eval(output));
        end
        

%         function output = apply(obj,input)
%             output = regexprep(input, ['\<',obj.name,'\>'], obj.value);  
%         end
        
%         function value = func(obj,~,~,~)
%             value = obj.value;
%         end

        
    end

end