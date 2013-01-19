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
            
            [idx, found] = obj.locate(name);
            

            
            if found;
                if any(obj.kernels(idx) == kern);
                    error('MatrixKernelRegistry:add', 'An instance of this kernel already exists in the registry.');
                end
                
%                 if ~isempty(kern.value);
%                     error('MatrixKernelRegistry:add','The matrix kernel being added already exists and has a non-empty Matrix.');
%                 end
                kern.matrix = obj.kernels(idx).matrix;
            end
            
            if isempty(obj.kernels)
                obj.kernels = kern;
            else
                obj.kernels(end+1) = kern;
            end
        end
        
        function K = assemble(obj,name,varargin)
            
           idx = obj.locate(name);
           
            for i = 1:length(idx);
                obj.kernels(idx(i)).assemble();
            end
           
           K = obj.kernels(idx(1)).matrix.init(); 
           
        end
        
        
        
    end
    
end

