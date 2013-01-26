classdef ConstantVectorRegistry < mFEM.registry.base.Registry
    %CONSTANTVECTORREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = protected)
        kernels = mFEM.kernels.ConstantVector.empty();
        options = struct('add', false, 'boundary', [], 'component', [], 'subdomain', []);
        mesh;
    end
    
    methods (Access = public)
        function obj = ConstantVectorRegistry(mesh, varargin)
            obj = obj@mFEM.registry.base.Registry(varargin{:});
            obj.mesh = mesh;
        end
        
         function value = get(obj, name)
                kern = obj.find(name);
                value = kern.value;
         end
    end 
      
    methods (Access = protected) 
        function knew = addKernel(obj, name, input, opt)

            if opt.add;
                % Create a new kernel
                knew = mFEM.kernels.ConstantVector(name, obj.mesh, input, ...
                    'Boundary', opt.boundary, 'Component', opt.component,...
                    'Subdomain', opt.subdomain);
                
                [idx, found] = obj.find(name, '-add');
                if found;
                    kold = obj.kernels(idx);                    
                    if kold.value.m ~= knew.value.m
                        error('ConstantVectorRegistry:addKernel', 'The vectors must be the same length');
                    end
                    kold.value.add(knew.value.init()); 
                    knew = kold;
                else
                    %k0 = mFEM.kernels.ConstantVector(name, input);
                    obj.kernels(idx) = knew;
                end
                           
            else
                idx = obj.find(name,'-add');      
                knew = mFEM.kernels.ConstantVector(name, obj.mesh, input, ...
                    'Boundary', opt.boundary, 'Component', opt.component,...
                    'Subdomain', opt.subdomain);
                obj.kernels(idx) = knew;
            end
        end
    end
end

