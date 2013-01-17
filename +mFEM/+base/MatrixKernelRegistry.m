classdef MatrixKernelRegistry < handle
    %CONSTANTKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        constants = mFEM.base.ConstantKernel.empty();
    end
    
    methods %(Access = Public)
        function obj = ConstantKernelRegistry()
        end
        
        function add(obj,name,value)
            obj.constants(end+1) = mFEM.base.ConstantKernel(name,value);
        end
        
        
    end
    
end

