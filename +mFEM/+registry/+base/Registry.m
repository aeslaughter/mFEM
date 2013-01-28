classdef Registry < handle
    properties
        reserved = {}% reserved variables, not available for constants
            %{'N','B','L','x','t','xi','eta','zeta','elem','Ke','grad'};
        core_options = struct('disablewarnings', false, 'allowduplicates', false);
    end
    
    properties (Abstract, GetAccess = public, SetAccess = protected)
        options;
        kernels;
    end
    
    methods (Abstract, Access = protected)
        kern = addKernel(obj, name, input, opt)
    end

   methods (Abstract)
        output = get(obj, name);
   end

    methods
               
        function obj = Registry(varargin)
            obj.core_options = gatherUserOptions(obj.core_options, varargin{:},'GatherUserOptions',{'-disablewarn'});
        end

        function str = apply(obj, str, varargin)        
            for i = 1:length(obj.kernels)   
                str = obj.kernels(i).apply(str, varargin{:});
            end                      
        end
        
        function kern = add(obj, varargin)
            [opt, input] = gatherUserOptions(obj.options, varargin{:});
            
            % Location of last input flag
            n = length(input) - 1;
            
            % Loop through each name and store in the const property
            for i = 1:2:n;
                name = input{i};
                value = input{i+1};
                kern = obj.addKernel(name, value, opt);
            end  
        end
        
        function [idx, found] = find(obj, name, varargin)

            %opt = obj.reg_options;
            opt.index = false;
            opt.add  = false;
            opt = gatherUserOptions(opt, varargin{:});
            
            % Locate the kernel
            idx = [];
            found = false;
            for i = 1:length(obj.kernels);
                if strcmp(name, obj.kernels(i).name)
                    idx = [idx,i];
                    found = true;
                end
            end
            
            % Found and adding a new kernel
            if opt.add && found;
%                 if ~obj.core_options.disablewarnings && ~obj.core_options.allowduplicates;
%                     warning('The variable %s was previously defined, the new value will replace the existing.', name);
% %                 end
%                 if ~obj.core_options.allowduplicates;
%                     return;
%                 end
                
            % Not found and adding a new kernel    
            elseif opt.add && ~found;
                idx = length(obj.kernels) + 1;
                
            % Found, not adding, and returning the actual kerenls    
            elseif ~opt.add && ~opt.index && found;
                idx = obj.kernels(idx); 
                
            % Found, not adding, and returning the indices (this is
            % default)
%             elseif opt.index && found;
                
            % Not found, not adding    
%             elseif ~found;
%                 error('Registry:find', 'The variable %s was not found.', name); 
            end
        end
    end          
end

