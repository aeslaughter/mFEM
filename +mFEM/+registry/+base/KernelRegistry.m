classdef KernelRegistry < handle
    properties
        reserved = ... % reserved variables, not available for constants
            {'N','B','L','x','t','xi','eta','zeta','elem','Ke','grad'};
        options = struct(...
            'type','constant','disablewarnings', false);
        kernels;
    end
    
    methods (Access = protected)
        function str = apply(obj, kern)
            str = kern.value;

            for i = 1:length(obj.kernels)   
                expr = obj.kernels(i).name;
                repstr = obj.kernels(i).value;
                str = regexprep(str,['\<',expr,'\>'],repstr);
            end                      
        end
    end
    
    methods
        function obj = KernelRegistry(varargin)

            
            obj.options = gather_user_options(obj.options,varargin{:});
   
            switch lower(obj.options.type);
                case {'c','const','constant'};
                    obj.kernels = mFEM.kernels.base.ConstantKernel.empty();
                case {'f','func','function'}
                    obj.kernels = mFEM.kernels.base.FunctionKernel.empty();
                case {'m','mat','matrix'}
                    obj.kernels = mFEM.kernels.base.MatrixKernel.empty();
                otherwise
                    error('KernelRegistry:add_kernel', 'The type %s was not recoqnized.', obj.options.type);
            end
            
            obj.options.type = class(obj.kernels);

        end
        
        function kern = add(obj,varargin)
            
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
            if any(strcmp(name,obj.reserved));
                error('KernelRegistry:test_name', 'The name %s is reserved and may not be used for naming variables of any type.', name);
            end
        end
  

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
      
    methods (Access = protected)
        function kern = add_kernel(obj,name,value)
            
            obj.test_name(name);
            [idx, found] = obj.locate(name);
   
            kern = feval(obj.options.type, name, value);

            if found && ~obj.options.disablewarnings;
                warning('The value %s was previously defined, the new value will replace the existing.', kern.name);
            end

            if isa(kern, 'mFEM.kernels.base.ConstantKernel');
                str = obj.apply(kern);
                kern.value = num2str(eval(str));
            end
            
            obj.kernels(idx) = kern;
        end  
    end
end

