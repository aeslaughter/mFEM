classdef KernelRegistry < handle
    properties
        reserved = ... % reserved variables, not available for constants
            {'N','B','L','x','t','xi','eta','zeta','elem','Ke','grad'};
        options = struct(...
            'disablewarnings', false);
    end
    
    properties (Abstract)
        kernels;
    end
    
    methods (Abstract, Access = protected)
        kern = add_kernel(obj,name,value,varargin)
    end
    
    methods
        function obj = KernelRegistry(varargin)
            obj.options = gather_user_options(obj.options,varargin{:});
        end
        
        function kern = add(obj,varargin)
            
            % special case
            if length(varargin) == 3;
                kern = obj.add_kernel(varargin{:});
                return;
            end
                   
            % Location of last input flag
            n = nargin - 2;
            
            % Loop through each name and store in the const property
            for i = 1:2:n;
                name = varargin{i};
                value = varargin{i+1};
                kern = obj.add_kernel(name, value);
            end  
            
        end
        function test_name(obj,name)
            
            %if isempty(obj); return; end
            if any(strcmp(name,obj.reserved));
                error('KernelRegistry:test_name', 'The name %s is reserved and may not be used for naming variables of any type.', name);
            end
        end
    end
        
    methods
        function [idx,found] = locate(obj, name)
            
            idx = [];
            found = false;
            for i = 1:length(obj.kernels);
                if strcmp(name, obj.kernels(i).name)
                    idx = i;
                    found = true;
                    return;
                end
            end
            
            if isempty(idx);
                idx = length(obj.kernels) + 1;
            end
            
        end
  
    end  
end

