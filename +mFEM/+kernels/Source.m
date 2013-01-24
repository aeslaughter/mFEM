classdef Source < mFEM.kernels.base.MatrixKernel
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
        function obj = Source(mesh,varargin)
            obj = obj@mFEM.kernels.base.MatrixKernel(mesh,'Source',varargin{:},'Type','vector');

            opt.b = 1;
            opt.function = false;
            
            [opt, ~] = gatherUserOptions(opt,varargin{:});
            
            if isa(opt.b,'function_handle') || (opt.function && ischar(opt.b));
                obj.B = mFEM.kernels.Func('b', opt.b);
            else
                obj.B = mFEM.kernels.Constant('b', opt.b);
            end
        end
        
        function value = eval(obj,elem,qp,varargin)
            N = @(qp) elem.shape(qp);
            x = elem.getPosition(qp);
            b = obj.B.eval(elem,x,varargin{:});
            value = N(qp)'*b;
        end
    end
    
end
