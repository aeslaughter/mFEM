classdef FunctionKernelRegistry < mFEM.registry.base.KernelRegistry
    %FUNCTIONKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        kernels = mFEM.kernels.base.FunctionKernel.empty();
        constReg;
        options = struct(...
            'disablewarnings', false, 'constants', []);
    end
    
    methods %(Access = Public)
        function obj = FunctionKernelRegistry(varargin)
            obj = obj@mFEM.registry.base.KernelRegistry(varargin{:});
            
            obj.options = gather_user_options(obj.options,varargin{:});

            if isa(obj.options.constants, 'mFEM.registry.ConstantKernelRegistry');
               obj.constReg = obj.options.constants;
            end
        end
    
        function kern = add(obj,name,input,varargin)
           
            opt.constants = obj.constReg;
            opt = gather_user_options(opt,varargin{:});
            kern = obj.add_kernel(name, input, 'constants', opt.constants);


        end
    end
end

