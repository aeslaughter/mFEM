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
        function obj = Diffusion(mesh, varargin)
            obj = obj@mFEM.kernels.base.MatrixKernel(mesh, 'Diffusion', 'Type', 'Matrix');

            opt.d = 1;
            opt.function = false;
            opt = gatherUserOptions(opt,varargin{:});
            
            if isa(opt.d,'function_handle') || (opt.function && ischar(opt.d));
                obj.D = mFEM.kernels.Func('D',opt.d);
            else
                obj.D = mFEM.kernels.Constant('D',opt.d);
            end
        end
        
        function value = eval(obj,elem,qp,varargin)
            B = @(qp) elem.shapeDeriv(qp);
            x = elem.getPosition(qp);
            D = obj.D.eval(elem,x,varargin{:});
            value = B(qp)'*D*B(qp);
        end
    end
    
end

