classdef FuncRegistry < mFEM.registry.base.Registry
    %FUNCREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = protected)
        kernels = mFEM.kernels.Func.empty();
        options = struct('constantregistry', mFEM.registry.ConstantRegistry.empty());
    end
    
    methods %(Access = Public)
        function obj = FuncRegistry(varargin)
            obj = obj@mFEM.registry.base.Registry(varargin{:});

            obj.options = gatherUserOptions(obj.options, varargin{:}, 'GatherUserOptoins', {'-disablewarn'});
            
        end
          
     function value = get(obj, name, varargin)
         
         kern = obj.find(name);
         
         if nargin == 2;
            if isa(obj.input, 'function_handle');
                value = kern.input;
            else
                value = str2func(kern.input);
            end
         elseif nargin == 4 || nargin == 5;
            value = kern.eval(varargin{:});
         else
             error('FuncRegistry:get:InvalidInput', 'Invalid number of inputs.');
         end
     end
     
    end
    
    methods (Access = protected)
        function kern = addKernel(obj, name, input, opt)
            
            idx = obj.find(name, '-add');
            kern = mFEM.kernels.Func(name, input,...
                'ConstantRegistry', opt.constantregistry);
         
            obj.kernels(idx) = kern;

       end   
        
        
    end
end

