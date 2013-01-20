classdef MatrixKernel < mFEM.kernels.base.Kernel ...
                      & matlab.mixin.Heterogeneous
    %MATRIXKERNEL Abstract class for defining finite element matrices

    properties
       mesh;
       options = struct(...
           'boundary', [], 'subdomain', []);
    end
    
    properties (SetAccess = {?mFEM.registry.MatrixKernelRegistry})
        matrix;
    end
   
    methods 
        function obj = MatrixKernel(mesh, name, varargin)
            obj = obj@mFEM.kernels.base.Kernel(name);

            obj.options = gather_user_options(obj.options, varargin{:});
            obj.mesh = mesh;
            obj.matrix = mFEM.Matrix(mesh);
        end

        function K = get(obj)
            K = obj.matrix.init();
        end
        
        function varargout = assemble(obj, varargin)
                
                opt.zero = false;
                opt = gather_user_options(opt, varargin{:});
                

               for e = 1:obj.mesh.n_elements;

                    elem = obj.mesh.element(e);
                    Ke = zeros(elem.n_dof);

                    % Loop over the quadrature points
                    for i = 1:length(elem.qp);
                        Ke = Ke + elem.W(i)*obj.eval(elem,elem.qp{i})*elem.detJ(elem.qp{i});
                    end

                    dof = elem.get_dof();
                    obj.matrix.add_matrix(Ke, dof); 
               end

               if nargout == 1;
                    varargout{1} = obj.matrix.init(); 
                    if opt.zero;
                        obj.matrix.zero();
                    end
               end
        end       
    end

end