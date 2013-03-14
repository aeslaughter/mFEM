classdef MatrixKernelRegistry < mFEM.base.Registry
    %MATRIXKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
   
    properties (GetAccess = public, SetAccess = protected)
        kernels;% = mFEM.kernels.base.MatrixKernel.empty();
        options = struct('Mesh', mFEM.FEmesh.empty());
        mesh;
    end
    
    methods %(Access = Public)
        function obj = MatrixKernelRegistry(mesh, varargin)
            obj = obj@mFEM.base.Registry(varargin{:},'-AllowDuplicates');
            obj.mesh = mesh;
        end 
        
%         function apply(~, varargin)
%            error('Not implmented'); 
%         end
        
        function value = get(obj, name)      
            kern = obj.find(name);
            value = kern.value;
         end
     
        function K = assemble(obj,name,varargin)
            
            opt.zero = false;
            opt.tag = {};
            opt.component = [];
            opt = gatherUserOptions(opt, varargin{:});   
                
            kern = obj.find(name);
           
            K = kern(1).assemble(varargin{:},'-parallel');
            for i = 2:length(kern);
                K = K + kern(i).assemble(varargin{:},'-parallel');
            end
           
            if opt.zero;
                for i = 1:length(kern);
                    kern(i).value.zero();
                end
            end           
           
        end
    end
    
    methods (Access = protected)
        function kern = addKernel(obj, name, input, varargin)
            
            obj.options.mesh = obj.mesh;
            opt = gatherUserOptions(obj.options,varargin{:});
            
            kern = input;
            kern.name = name;          
            
            [idx, found] = obj.find(name,'-add');

            if found;
                if any(obj.kernels(idx) == kern);
                    error('MatrixKernelRegistry:addKernel', 'An instance of this kernel already exists in the registry.');
                end
                
%                 if ~isempty(kern.value);
%                     error('MatrixKernelRegistry:add','The matrix kernel being added already exists and has a non-empty Matrix.');
%                 end
                kern.value = obj.kernels(idx).value;
            end
            
            if isempty(obj.kernels)
                obj.kernels = kern;
            else
                obj.kernels(end+1) = kern;
            end
        end
        
    end
    
end

