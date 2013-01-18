classdef Diffusion < mFEM.kernels.base.MatrixKernel
    %DIFFUSION A kernel for the diffusion equation.
    %   
    %
    % \nabla(w)*D*\nabla(u)
    %   
    % B'*D*B
    
    properties
       D;% = mFEM.base.ConstantKernel.empty(); 
        
    end
    
    methods
        function obj = Diffusion(varargin)
            
            
            obj = obj@mFEM.kernels.base.MatrixKernel('Diffusion');

            opt.d = 1;
            opt.function = false;
            
            
            opt = gather_user_options(opt,varargin{:});
            
            if isa(opt.d,'function_handle') || (opt.function && ischar(opt.d));
                obj.D = mFEM.kernels.base.FunctionKernel('D',opt.d);
            else
                obj.D = mFEM.kernels.base.ConstantKernel('D',opt.d);
            end
        end
        
        function value = eval(obj,elem,qp,varargin)
            B = @(qp) elem.shape_deriv(qp);
            x = elem.get_position(qp);
            D = obj.D.eval(elem,x,varargin{:});
            value = B(qp)'*D*B(qp);
        end
    end
    
end

