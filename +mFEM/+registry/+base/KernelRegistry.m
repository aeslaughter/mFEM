classdef KernelRegistry < handle
    properties
        reserved = ... % reserved variables, not available for constants
            {'N','B','L','x','t','xi','eta','zeta','elem','Ke','grad'};
        options = struct('disablewarnings', false, 'allowduplicates',false);
    end
    
    properties (Abstract)
        kernels;
    end
    
    methods (Abstract)
        kern = add(obj, varargin)
        varargout = apply(obj, varargin)
    end

    methods
        function obj = KernelRegistry(varargin)
            obj.options = gather_user_options(obj.options,varargin{:},{'-disablewarn'});
        end
        
%         function apply(obj, kern)
%             for i = 1:length(obj.kernels)   
%                 obj.kernels(i).apply(kern);
%             end                      
%         end

        function test_name(obj,name)
            if any(strcmp(name,obj.reserved));
                error('KernelRegistry:test_name', 'The name %s is reserved and may not be used for naming variables of any type.', name);
            end
        end
  
        function [idx, found] = locate(obj, name, varargin)

            opt = obj.options;
            opt.index = false;
            opt = gather_user_options(opt, varargin{:});
            
            
            idx = [];
            found = false;
            for i = 1:length(obj.kernels);
                if strcmp(name, obj.kernels(i).name)
                    idx = [idx,i];
                    found = true;
                    if ~obj.options.disablewarnings && ~obj.options.allowduplicates;
                        warning('The value %s was previously defined, the new value will replace the existing.', name);
                    end
                    if ~obj.options.allowduplicates;
                        return;
                    end
                end
            end
            
            if ~opt.index;
                if ~isempty(idx);
                    idx = obj.kernels(idx);
                end
                
            elseif isempty(idx);
                idx = length(obj.kernels) + 1;
            end
        end
    end          
end

