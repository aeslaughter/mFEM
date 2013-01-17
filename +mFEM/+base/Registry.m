classdef Registry < mFEM.base.ConstantKernelRegistry ...
                  & mFEM.base.MatrixKernelRegistry
    properties
%         reserved = ... % reserved variables, not available for constants
%             {'N','B','L','x','t','xi','eta','zeta','elem','Ke','grad'};
        

    end
    
    methods
        function obj = Registry()
 
        end
        
        
        function add_constant(obj, varargin)
            obj.add('constant',varargin{:});
        end
        
        function add_matrix(obj, varargin)
           % obj.add@mFEM.base.MatrixKernelRegistry(varargin{:});
            
        end

    end
    
    methods %(Access = private)
        
        function add(obj,type,varargin)
            switch lower(type)
                case 'constant'
                    obj.add@mFEM.base.ConstantKernelRegistry(varargin{:});
            end
    
        end
        
    end
    
end

