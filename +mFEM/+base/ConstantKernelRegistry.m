classdef ConstantKernelRegistry < mFEM.base.KernelRegistry
    %CONSTANTKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        constants = mFEM.base.ConstantKernel.empty();
    end
    
    methods %(Access = Public)
        function obj = ConstantKernelRegistry()
        end
        
        function add(obj,name,value,varargin)
            
            opt.add = false;
            opt = gather_user_options(opt,varargin{:});
            
            obj.test_name(name);
            [idx,found] = obj.locate(name,obj.constants);
            
            new_const = mFEM.base.ConstantKernel(name,value);
            obj.apply(new_const);

            if found && opt.add;    
                obj.constants(idx).merge(new_const);
            else
                obj.constants(idx) = new_const;  
            end
        
        end
        
        function apply(obj, kernel)
            
            for i = 1:length(obj.constants)
                kernel.input = obj.constants(i).apply(kernel.input);
            end
            
            kernel.input = num2str(eval(kernel.input));
            
        end
        
        
    end
    
end

