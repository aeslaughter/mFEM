classdef ConstantKernelRegistry < mFEM.registry.base.KernelRegistry
    %CONSTANTKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        kernels = mFEM.kernels.base.ConstantKernel.empty();
    end
    
    methods %(Access = Public)
        function obj = ConstantKernelRegistry(varargin)
            obj = obj@mFEM.registry.base.KernelRegistry(varargin{:});
        end

        function str = apply(obj, str)
            for i = 1:length(obj.kernels)
                str = obj.kernels(i).apply(str);
            end                      
        end
        
    end
    
    methods (Access = protected)
        function kern = add_kernel(obj,name,value,varargin)
                   
            opt.add = false;
            opt = gather_user_options(opt,varargin{:});

            obj.test_name(name);
            [idx, found] = obj.locate(name);
            
            kern = mFEM.kernels.base.ConstantKernel(name,value);

            if found && opt.add;   
                if ~obj.options.disablewarnings;
                    warning('The constant %s was previously defined, the new value is being added to the existing.', name);
                end
                kern.merge(obj.kernels(idx));
                
            elseif found
                if ~obj.options.disablewarnings;
                    warning('The constant %s was previously defined, the new value will replace the existing.', name);
                end
            end

            str = obj.apply(kern.value);
            kern.value = num2str(eval(str));
            obj.kernels(idx) = kern;
            

        end
        
    end
    
end

