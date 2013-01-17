classdef Kernel < handle;%< matlab.mixin.Heterogeneous
    properties(Access = public);
        name;
        %t;
    end
      
    methods 
        function obj = Kernel(name,varargin)
            obj.name = name;
        end
    end
    
    methods (Abstract)
        value = eval(obj,elem,qp,t);
        %str = apply(obj,str);
    end
end

