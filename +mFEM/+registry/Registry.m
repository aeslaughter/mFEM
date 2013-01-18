classdef Registry %< handle 
    properties
        constants = mFEM.registry.base.ConstantKernelRegistry();
        functions = mFEM.registry.base.FunctionKernelRegistry();
    end
    
    methods
        function obj = Registry()
 
        end
            
%         function delete(obj)
%            delete(obj.constants);
%            delete(obj.functions);
%         end
        
        function add_constant(obj, varargin)
            obj.constants.add(varargin{:});
        end
        
        function add_function(obj, varargin)
            kern = obj.functions.add(varargin{:});
            
            if ischar(kern.input);
                kern.input = obj.constants.apply(kern.input);
            end
        end
        
        function add_matrix(obj, varargin)
           % obj.add@mFEM.base.MatrixKernelRegistry(varargin{:});
            
        end

    end
    
%     methods %(Access = private)
%         
%         function add(obj,type,varargin)
%             switch lower(type)
%                 case 'constant'
%                     obj.add@mFEM.kernels.base.ConstantKernelRegistry(varargin{:});
%             end
%     
%         end
%         
%     end
    
end

