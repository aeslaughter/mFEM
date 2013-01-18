classdef MatrixKernel < mFEM.kernels.base.Kernel;

    properties
       options = struct(...
           'boundary',[],'subdomain',[]);
       mesh;
       matrix;
    end
    
    
    methods 
        function obj = MatrixKernel(mesh, name, varargin)
            obj = obj@mFEM.kernels.base.Kernel(name);
            obj.options = gather_user_options(obj.options, varargin{:});
            obj.mesh = mesh;
            obj.value = mFEM.Matrix(mesh);
        end
        
        function K = assemble(obj)
            K = obj.value.init();
        end
        
        
        function integrate(obj,e)
            % Initialize the stiffness matrix (K) and the force vector (f), for
            % larger systems K should be sparse.
            
            elem = obj.mesh.element(e);
            
            Ke = zeros(elem.n_dof);

            % Loop over the quadrature points in the two dimensions to perform the
            % numeric integration
            for i = 1:length(elem.qp);
                % Build stiffness matrix
                Ke = Ke + elem.W(i)*obj.eval(elem,elem.qp{i})*elem.detJ(elem.qp{i});
            end
            
            dof = elem.get_dof();
            obj.matrix.add_matrix(Ke, dof);
            
        end
        
    end

end