classdef Kernel < handle
    properties(Access = public);
        name;
        value;
    end
      
    methods 
        function obj = Kernel(name)
            obj.name = name;
        end    
    end
    
    methods (Abstract)
        value = eval(obj,elem,qp,t);

    end
    
%     methods (Access = private)
%         function str = apply(obj,str,varargin)
%             % Apply OBJ's value to STR
%             if ischar(obj.value);
%                 str = regexprep(obj.value, ['\<',obj.name,'\>'], str); 
%             end
%         end
%     end
end

