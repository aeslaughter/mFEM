classdef FunctionKernelRegistry < mFEM.registry.base.KernelRegistry
    %FUNCTIONKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        kernels = mFEM.kernels.base.FunctionKernel.empty();
    end
    
    methods %(Access = Public)
        function obj = FunctionKernelRegistry(varargin)
            obj = obj@mFEM.registry.base.KernelRegistry(varargin{:});
        end
        
%         function apply(obj, kernel)
%             kernel.input = func2str(kernel.input);
%             for i = 1:length(obj.kernels)
%                 kernel.input = obj.kernels(i).apply(kernel.input);
%             end
%             
%             %kernel.value = num2str(eval(kernel.value));
%             
%         end
        
    end
    
    methods (Access = protected)
        function kern = add_kernel(obj,name,value,varargin)
                   
            opt.add = false;
            opt = gather_user_options(opt,varargin{:});
           
            obj.test_name(name);
            [idx,found] = obj.locate(name);
            
            kern = mFEM.kernels.base.FunctionKernel(name,value);

            if found && opt.add;   
                error('The function %s was previously defined, adding functions is not yet supported.', name);
                
            elseif found
                if ~obj.options.disablewarnings;
                    warning('The function %s was previously defined, the new value will replace the existing.', name);
                end
            end
            
            obj.kernels(idx) = kern;
        end
        
    end
    
end

