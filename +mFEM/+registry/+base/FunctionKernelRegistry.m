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
    end
    

end

