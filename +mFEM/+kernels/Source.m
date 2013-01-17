classdef Source < mFEM.base.Kernel
    %SOURCE A kernel for the diffusion equation.
    %   
    %
    % w*b
    %   
    % N'*b
    
    properties
       B;
        
    end
    
    methods
        function obj = Source(varargin)
            obj = obj@mFEM.base.Kernel('Source');

            opt.b = 1;
            opt.function = false;
            
            opt = gather_user_options(opt,varargin{:});
            
            if isa(opt.b,'function_handle') || (opt.function && ischar(opt.b));
                obj.B = mFEM.base.FunctionKernel('b',opt.b);
            else
                obj.B = mFEM.base.ConstantKernel('b',opt.b);
            end
        end
        
        function value = eval(obj,elem,qp,varargin)
            N = @(qp) elem.shape(qp);
            x = elem.get_position(qp);
            b = obj.B.eval(elem,x,varargin{:});
            value = N(qp)'*b;
        end
    end
    
end
