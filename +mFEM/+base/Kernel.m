classdef Kernel %< matlab.mixin.Heterogeneous
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
        function obj = Kernel(system, name, input)
            obj.system = system;
            obj.name = name;
            obj.input = input;
            obj.value = obj.init();
        end
        
        function output = init(obj)
            
            output = obj.input;
            for i = 1:length(obj.system.kernels);
                output = obj.system.kernels{i}.apply(output);
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

