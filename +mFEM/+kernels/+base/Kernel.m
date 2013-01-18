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

end

