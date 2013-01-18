classdef MatrixKernelRegistry < mFEM.registry.base.KernelRegistry
    %CONSTANTKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
   
    properties
        kernels = mFEM.kernels.base.ConstantKernel.empty();
    end
    
    methods %(Access = Public)
        function obj = MatrixKernelRegistry(varargin)
            obj = obj@mFEM.registry.base.KernelRegistry(varargin{:});
        end  
    end
    
end

