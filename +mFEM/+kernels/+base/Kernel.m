classdef Kernel < handle
    properties(Access = public);
        name;
        input;
        %t;
    end
      
    methods 
        function obj = Kernel(name,input)
            obj.name = name;
            obj.input = input;
        end
        
%         function str = apply(obj,str)
%             str = regexprep(str, ['\<',obj.name,'\>'], obj.value);  
%         end
        
    end
    
    methods (Abstract)
        value = eval(obj,elem,qp,t);
        %str = apply(obj,varargin);
    end
end

