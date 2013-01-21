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
        
        function kern = add(obj,varargin)
            
            % Location of last input flag
            n = nargin - 2;
            
            % Loop through each name and store in the const property
            for i = 1:2:n;
                name = varargin{i};
                value = varargin{i+1};
                kern = obj.add_kernel(name, value);
            end  
        end

        function apply(obj, kern)
            for i = 1:length(obj.kernels)   
                obj.kernels(i).apply(kern);
            end                      
        end
    end 
      
    methods (Access = protected)
        function kern = add_kernel(obj,name,value,varargin)
            
            obj.test_name(name);
            idx = obj.locate(name);
   
            kern = mFEM.kernels.base.ConstantKernel(name, value, varargin{:});

            obj.apply(kern);

            obj.kernels(idx) = kern;
        end
    end
end

