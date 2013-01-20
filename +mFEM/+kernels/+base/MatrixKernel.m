classdef MatrixKernel < mFEM.kernels.base.Kernel ...
                      & matlab.mixin.Heterogeneous
    %MATRIXKERNEL Abstract class for defining finite element matrices

    properties
       mesh;
       options = struct(...
           'boundary', [], 'subdomain', [],'type', 'matrix');
    end
    
    properties (SetAccess = {?mFEM.registry.MatrixKernelRegistry})
        matrix;
    end
   
    methods 
        function obj = MatrixKernel(mesh, name, varargin)
            obj = obj@mFEM.kernels.base.Kernel(name);

            obj.options = gather_user_options(obj.options, varargin{:});
            obj.mesh = mesh;
            
            if any(strcmpi(obj.options.type,{'matrix','mat','m'}));
                obj.matrix = mFEM.Matrix(mesh);
                obj.options.type = 'matrix';
            elseif any(strcmpi(obj.options.type,{'vector','vec','v'}));
                obj.matrix = mFEM.Vector(mesh);
                obj.options.type = 'vector';
            else
                error('MatrixKernel:MatrixKernel', 'Unknown type %s.', obj.options.type);
            end
        end

        function K = get(obj)
            K = obj.matrix.init();
        end
        
        function varargout = assemble(obj, varargin)
                
                opt.zero = false;
                opt = gather_user_options(opt, varargin{:});
                

                
                
                
               for e = 1:obj.mesh.n_elements;

                    elem = obj.mesh.element(e);
                    
                    if strcmpi(obj.options.type,'matrix');
                          Ke = zeros(elem.n_dof);
                    else
                          Ke = zeros(elem.n_dof,1);         
                    end

                    % Loop over the quadrature points
                    for i = 1:length(elem.qp);
                        Ke = Ke + elem.W(i)*obj.eval(elem,elem.qp{i})*elem.detJ(elem.qp{i});
                    end

                    dof = elem.get_dof();
                    obj.matrix.add(Ke, dof); 
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