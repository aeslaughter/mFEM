classdef Registry < handle
    properties
        reserved = {}% reserved variables, not available for constants
            %{'N','B','L','x','t','xi','eta','zeta','elem','Ke','grad'};
        options = struct('disablewarnings', false, 'allowduplicates',false);
    end
    
%     properties (Access = protected)
%         kernels;
%     end
    
%    methods (Abstract)
%          output = findKernel(obj, varargin)
%         kern = add(obj, varargin)
%         varargout = apply(obj, varargin)
%    end

    methods
        function obj = Registry(varargin)
            obj.options = gather_user_options(obj.options,varargin{:},{'-disablewarn'});
        end

        function testName(obj,name)
            if any(strcmp(name,obj.reserved));
                error('KernelRegistry:test_name', 'The name %s is reserved and may not be used for naming variables of any type.', name);
            end
        end
  
        function [idx, found] = findKernel(obj, name, kernels, varargin)

            opt = obj.options;
            opt.index = false;
            opt = gather_user_options(opt, varargin{:});

            idx = [];
            found = false;
            for i = 1:length(kernels);
                if strcmp(name, kernels(i).name)
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
                    idx = kernels(idx);
                end
                
            elseif isempty(idx);
                idx = length(kernels) + 1;
            end
        end
    end          
end

