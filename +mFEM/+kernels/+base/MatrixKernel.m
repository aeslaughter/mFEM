classdef MatrixKernel < mFEM.kernels.base.Kernel                     
    %MATRIXKERNEL Abstract class for defining finite element matrices

    properties
       options = struct(...
           'boundary', [], 'subdomain', [], 'component', [], 'type', 'matrix');
       reserved = {'N','B','Ke'};
    end
    
    properties (SetAccess = {?mFEM.registry.base.Registry})
        value;
    end
   
    properties (Access = protected)
    	mesh;
        direct = false;
    end
       
    methods 
        function obj = MatrixKernel(mesh, name, varargin)
            obj = obj@mFEM.kernels.base.Kernel(name);

            [obj.options, unknown] = gatherUserOptions(obj.options, varargin{:});
            obj.mesh = mesh;
            
            if any(strcmpi(obj.options.type,{'matrix','mat','m'}));
                obj.value = mFEM.Matrix(mesh);
                obj.options.type = 'matrix';
            elseif any(strcmpi(obj.options.type,{'vector','vec','v'}));
                obj.value = mFEM.Vector(mesh);
                obj.options.type = 'vector';
            else
                error('MatrixKernel:MatrixKernel', 'Unknown type %s.', obj.options.type);
            end
        end

        function str = apply(obj, str, elem, qp, t)
            error('MatrixKernel:apply', 'Not implemented');
%             if ~ischar(str);
%                 error('Func:apply', 'The input (str) must be a character string');
%             end    
%             
%             % Apply OBJ's value to KERN
%             expr = ['\<',obj.name,'\>'];
%             repstr = obj.eval(elem, qp, t);
%             
%             str = regexprep(str, expr, repstr); 
        end
        
        function K = get(obj)
            K = obj.value.init();
        end
        
        function varargout = assemble(obj, varargin)
                
                opt.zero = false;
                opt.boundary = obj.options.boundary;
                opt.subdomain = obj.options.subdomain;
                opt.component = obj.options.component;
                opt = gatherUserOptions(opt, varargin{:});

                elem = obj.mesh.getElements('boundary', opt.boundary, ...
                                             'subdomain', opt.subdomain);
                                            %'component', opt.component);
                
               for i = 1:length(elem);
               
                    if isempty(opt.boundary);
                        Ke = obj.evaluateElement(elem(i));
                    else
                        Ke = obj.evaluateSide(elem(i), opt.boundary);
                    end

                    dof = elem(i).getDof();
                    obj.value.add(Ke, dof); 
               end

               if nargout == 1;
                    varargout{1} = obj.value.init(); 
                    if opt.zero;
                        obj.value.zero();
                    end
               end
        end   
        
    end
    
    methods (Access = protected)
        function Ke = evaluateElement(obj, elem)
              
            if obj.direct
                error('Not yet supported');
                %Ke = obj.eval(elem,elem.qp{i})*elem.detJ(elem.qp{i});
                return;
            end
            
            if strcmpi(obj.options.type,'matrix');
                  Ke = zeros(elem.n_dof);
            else
                  Ke = zeros(elem.n_dof,1);         
            end

            % Loop over the quadrature points
            for i = 1:length(elem.qp);
                Ke = Ke + elem.W(i)*obj.eval(elem, elem.qp{i})*elem.detJ(elem.qp{i});
            end
        end
        
        function Ke = evaluateSide(obj, elem, id)
        
            if strcmpi(obj.options.type,'matrix');
                  Ke = zeros(elem.n_dof);
            else
                  Ke = zeros(elem.n_dof,1);         
            end
            
            if elem.local_n_dim == 1;
                for s = 1:elem.n_sides; 
                    if any(elem.side(s).boundary_id == id);
                        side = elem.buildSide(s);
                        dof = elem.getDof('Side', s, '-local');
                        Ke(dof) = Ke(dof) + obj.eval(side,[]);  
                        delete(side);
                    end
                end     
            else
                error('2D side not working yet');
            end

        end
    end

end