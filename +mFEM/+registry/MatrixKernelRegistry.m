classdef MatrixKernelRegistry < mFEM.registry.base.KernelRegistry
    %MATRIXKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
   
    properties
        kernels;% = mFEM.kernels.base.MatrixKernel.empty();
        mesh;
    end
    
    methods %(Access = Public)
        function obj = MatrixKernelRegistry(mesh,varargin)
            obj = obj@mFEM.registry.base.KernelRegistry(varargin{:},'-AllowDuplicates');
            obj.mesh = mesh;
        end 
        
        function kern = add(obj,name,input,varargin)
            opt.mesh = obj.mesh;
            opt = gather_user_options(opt,varargin{:});
            
            kern = input;
            kern.name = name;
            kern
            
            if isempty(obj.kernels)
                obj.kernels = kern;
            else
                obj.kernels(end+1) = kern;
            end
            
            
        end
        
    end
    
end

