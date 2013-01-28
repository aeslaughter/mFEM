classdef Func < mFEM.kernels.base.Kernel;
    properties %(Access = protected)
        input;  % function as input
        reserved = {'elem', 'x', 't'};
        const;
    end

    methods 
        function obj = Func(name, input, varargin)
            obj = obj@mFEM.kernels.base.Kernel(name);
            
            opt.constantregistry = mFEM.registry.ConstantRegistry.empty();
            opt.reserved = {};
            opt = gatherUserOptions(opt, varargin{:});
            
            obj.const = opt.constantregistry;
            obj.reserved = [obj.reserved, opt.reserved];
            obj.testName(name);
  
            if isa(input,'function_handle') 
                n = nargin(input);
                if n ~= 3;
                    error('Func:Func:InvalidHandle', 'The input function supplied accepts %d agruments, but it must accept three arguments: elem, x, and t.', n);
                end
            elseif ~ischar(input);
                error('Func:Func:InvalidInput', 'The input must be a function handle or a character string convertable to a handle, but a %s was given.', class(input));    
            end

            obj.input = input;
                    
        end        
        
        function fhandle = get(obj)
            if ishcar(obj.input)
                fhandle = str2func(obj.input);
            else
                fhandle = obj.input;
            end
        end
        
        function str = apply(obj, str, elem, x, t)
            
            if ~ischar(str);
                error('Func:apply', 'The input (str) must be a character string');
            end    
            
            % Apply OBJ's value to KERN
            expr = ['\<',obj.name,'\>'];
            repstr = mat2str(obj.eval(elem, x, t));
            
            str = regexprep(str, expr, repstr); 
        end
        
        function value = eval(obj,elem,x,varargin)
            t = [];
            if nargin == 4; 
                t = varargin{1};
            end
            
            if ischar(obj.input);
                if ~isempty(obj.const);
                    str = obj.const.apply(obj.input);
                else
                    str = obj.input;
                end
            
                fcn = str2func(['@(elem,x,t) ', str]);
            else
                fcn = obj.input;
            end

            value = feval(fcn, elem, x, t);
        end

        
    end

end