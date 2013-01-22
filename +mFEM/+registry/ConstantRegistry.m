classdef ConstantRegistry < mFEM.registry.base.Registry
    %CONSTANTKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        const = mFEM.kernels.Constant.empty();
    end
    
    methods %(Access = Public)
        function obj = ConstantRegistry(varargin)
            obj = obj@mFEM.registry.base.Registry(varargin{:});
        end
        
        function kern = addConstant(obj,varargin)
            
            % Location of last input flag
            n = nargin - 2;
            
            % Loop through each name and store in the const property
            for i = 1:2:n;
                name = varargin{i};
                value = varargin{i+1};
                kern = obj.addConstantKernel(name, value);
            end  
        end

        function str = applyConstant(obj, str)
            
            for i = 1:length(obj.const)   
                
                % Apply OBJ's value to str
                expr = ['\<',obj.const(i).name,'\>'];
                repstr = obj.const(i).value;

                if ~ischar(str);
                    error('ConstantKernel:applyConstant', 'The supplied function must be a character string');
                end 

                str = regexprep(str, expr, repstr); 

            end                      
        end

        function idx = findConstant(obj, name, varargin)
            
            idx = obj.findKernel(name, obj.const, varargin{:});

        end

    end 
      
    methods (Access = protected)
        
        function kern = addConstantKernel(obj,name,value,varargin)
            
            obj.testName(name);
            idx = obj.findConstant(name,'-index');
   
            kern = mFEM.kernels.Constant(name, value, varargin{:});

            kern.value = obj.applyConstant(kern.value);

            obj.const(idx) = kern;
        end
    end
end

