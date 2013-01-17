classdef KernelRegistry < handle
    properties
        reserved = ... % reserved variables, not available for constants
            {'N','B','L','x','t','xi','eta','zeta','elem','Ke','grad'};
        

    end
    
    methods
        function obj = KernelRegistry()

        end
        
        function test_name(obj,name)
            
            if any(strcmp(name,obj.reserved));
                error('KernelRegistry:test_name', 'The name %s is reserved and may not be used for naming variables of any type.', name);
            end
        end
    end
        
    methods (Static)
        function [idx,found] = locate(name,kernels)
            
            idx = [];
            found = false;
            for i = 1:length(kernels);
                if strcmp(name, kernels(i).name)
                    idx = i;
                    found = true;
                    return;
                end
            end
            
            if isempty(idx);
                idx = length(kernels) + 1;
            end
            
        end
        
        
        
    end
    
end

