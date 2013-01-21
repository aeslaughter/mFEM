classdef FunctionKernelRegistry < mFEM.registry.base.KernelRegistry
    %FUNCTIONKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        kernels = mFEM.kernels.base.FunctionKernel.empty();
        constReg;
    end
    
    methods %(Access = Public)
        function obj = FunctionKernelRegistry(varargin)
            obj = obj@mFEM.registry.base.KernelRegistry(varargin{:});
            
            obj.options.constants = [];  
            obj.options = gather_user_options(obj.options,varargin{:});

            if isa(obj.options.constants, 'mFEM.registry.ConstantKernelRegistry');
               obj.constReg = obj.options.constants;
            end
        end
        
        function output = apply(obj, kern, input, elem, x, t)
            
            output = input;
            for i = 1:length(obj.kernels)   
                output = obj.kernels(i).apply(kern, output, elem, x, t);
            end                      
        end
        
        function kern = add(obj,name,input,varargin)
           
            opt.constants = obj.constReg;
            opt = gather_user_options(opt,varargin{:});

            obj.test_name(name);
            idx = obj.locate(name);
   
            kern = mFEM.kernels.base.FunctionKernel(name, input, 'constants', opt.constants);

            obj.kernels(idx) = kern;

        end
    end
end

