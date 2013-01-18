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
        kern = add(obj,varargin)
    end

    methods
        function obj = KernelRegistry(varargin)
            obj.options = gather_user_options(obj.options,varargin{:},{'-disablewarn'});
        end
        
        function apply(obj, kern)
            for i = 1:length(obj.kernels)   
                obj.kernels(i).apply(kern);
            end                      
        end

        function test_name(obj,name)
            if any(strcmp(name,obj.reserved));
                error('KernelRegistry:test_name', 'The name %s is reserved and may not be used for naming variables of any type.', name);
            end
        end
  
        function idx = locate(obj, name)

            idx = [];
            for i = 1:length(obj.kernels);
                if strcmp(name, obj.kernels(i).name)
                    idx = [idx,i];
                    
                    if ~obj.options.disablewarnings && ~obj.options.allowduplicates;
                        warning('The value %s was previously defined, the new value will replace the existing.', name);
                    end
                    if ~obj.options.allowduplicates;
                        return;
                    end
                end
            end
            
            if isempty(idx);
                idx = length(obj.kernels) + 1;
            end
        end
    end    
      
%     methods (Access = protected)
%         function kern = add_kernel(obj,name,value,varargin)
%             
%             obj.test_name(name);
%             idx = obj.locate(name);
%    
%             kern = feval(obj.type, name, value, varargin{:});
% 
%             if isa(kern, 'mFEM.kernels.base.ConstantKernel');
%                 obj.apply(kern);
%             end
%             
%             obj.kernels(idx) = kern;
%         end  
%    end
end

