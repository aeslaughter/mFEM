classdef ConstantRegistry < mFEM.base.Registry
    %CONSTANTKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = protected)
        kernels = mFEM.kernels.Constant.empty();
        options = struct('add', false);
    end
    
    methods (Access = public)
        function obj = ConstantRegistry(varargin)
            obj = obj@mFEM.base.Registry(varargin{:});
        end
        
     function value = get(obj, name)
            kern = obj.find(name);
            value = kern.eval();
        end
        
    end 
      
    methods (Access = protected) 
        function k0 = addKernel(obj,name,value,opt)

            if opt.add;
                k1 = mFEM.kernels.Constant(name, value);
                k1.value = obj.apply(k1.value);
                
                [idx, found] = obj.find(name, '-add');
                if found;
                    k0 = obj.kernels(idx);
                    k0.value = [k0.value,'+',k1.value]; 
                else
                    k0 = mFEM.kernels.Constant(name, value);
                    obj.kernels(idx) = k0;
                end
                           
            else
                idx = obj.find(name,'-add');      
                k0 = mFEM.kernels.Constant(name, value);
                k0.value = obj.apply(k0.value);
                obj.kernels(idx) = k0;
            end
        end
    end
end

