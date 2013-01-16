classdef ConstantKernel < mFEM.base.Kernel;
%     properties(Access = public);
%         value;
% 
%     end

        
    methods 
        function obj = ConstantKernel(name, input)
            obj = obj@mFEM.base.Kernel(name, input);
            if isnumeric(obj.input);
                obj.input = num2str(obj.input);
            end
        end
        
        function output = init(obj,kernels)
            output = init@mFEM.base.Kernel(obj,kernels);
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