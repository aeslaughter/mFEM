classdef ConstantKernelRegistry < mFEM.base.KernelRegistry
    %CONSTANTKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        constants = mFEM.base.ConstantKernel.empty();
    end
    
    methods %(Access = Public)
        function obj = ConstantKernelRegistry(varargin)
            obj = obj@mFEM.base.KernelRegistry(varargin{:});
        end
        
        function add(obj,varargin)
            
            % special case
            if length(varargin) == 3;
                obj.add_private(varargin{:})
                return;
            end
                   
            % Location of last input flag
            n = nargin - 2;
            
            % Loop through each name and store in the const property
            for i = 1:2:n;
                name = varargin{i};
                value = varargin{i+1};
                obj.add_private(name,value);
            end  
            
        end
        
        function apply(obj, kernel)
            kernel.value = kernel.input;
            for i = 1:length(obj.constants)
                kernel.value = obj.constants(i).apply(kernel.value);
            end
            
            kernel.value = num2str(eval(kernel.value));
            
        end
        
    end
    
    methods (Access = private)
        function add_private(obj,name,value,varargin)
                   
            opt.add = false;
            opt = gather_user_options(opt,varargin{:});
           
            
            obj.test_name(name);
            [idx,found] = obj.locate(name,obj.constants);
            
            new_const = mFEM.base.ConstantKernel(name,value);

            if found && opt.add;   
                if ~obj.options.disablewarnings;
                    warning('The constant %s was previously defined, the new value is being added to the existing.', name);
                end
                new_const.merge(obj.constants(idx));
            elseif found
                if ~obj.options.disablewarnings;
                    warning('The constant %s was previously defined, the new value will replace the existing.', name);
                end
            end
            
            obj.constants(idx) = new_const;
            obj.apply(new_const);

        end
        
    end
    
end

