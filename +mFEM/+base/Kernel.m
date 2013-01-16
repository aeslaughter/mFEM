classdef Kernel < handle;%< matlab.mixin.Heterogeneous
    properties(Access = public);
        system;
        name;
        input;
        value;
       % type;
        
%         elem_;
%         x_;
%         t_;
    end

%     properties (Abstract)
%         value;
%     end
        
    methods 
        function obj = Kernel(name, input)
            obj.name = name;
            obj.input = input;
        end
        
        function output = init(obj,kernels)
            output = obj.input;
            for i = 1:length(kernels);
                output = kernels{i}.apply(output);
            end 
        end
        
        function output = apply(obj,input)
            output = regexprep(input, ['\<',obj.name,'\>'], obj.value);  
        end
        
        
    end
    
%     methods (Abstract)
%         output = apply(obj,input);
%     end
end

