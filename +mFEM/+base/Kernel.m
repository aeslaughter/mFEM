classdef Kernel < handle
    properties(Access = public);
        name;
    end
      
    properties (Abstract)%, Access = protected)
        reserved;
    end
    
    methods 
        function obj = Kernel(name)
            obj.name = name;
        end  
        
        function testName(obj,name)
            if any(strcmp(name,obj.reserved));
                error('Kernel:testName:InvalidName', 'The name %s is reserved and may not be used.', name);
            end
        end
    end
    
    methods (Abstract)
        value = eval(obj,elem,qp,t);
        str = apply(obj, str, varargin);
        value = get(obj);
    end

end

