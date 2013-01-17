classdef Registry < mFEM.base.ConstantKernelRegistry
    properties
%         reserved = ... % reserved variables, not available for constants
%             {'N','B','L','x','t','xi','eta','zeta','elem','Ke','grad'};
        

    end
    
    methods
        function obj = Registry()
            
            
        end
        
        
        function add_constant(obj, varargin)
        
            % special case
            if length(varargin) == 3;
                obj.add(varargin{:})
                return;
            end
                   
            % Location of last input flag
            n = nargin - 2;
            
            % Loop through each name and store in the const property
            for i = 1:2:n;
                name = varargin{i};
                value = varargin{i+1};
                obj.add(name,value);
            end 
        end

    end
    
end

