classdef MatrixKernel < mFEM.kernels.base.Kernel                     
    %MATRIXKERNEL Abstract class for defining finite element matrices

    properties
       t = 0;
       options = struct(...
           'tag', [], 'component', [], 'type', 'matrix');
       reserved = {'N','B','Ke','elem','qp','x','t','L'};

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

        function apply(~,varargin)
            error('MatrixKernel:apply', 'Not implemented.');
        end
        
        function K = get(obj,varargin)
            
            opt.init = false;
            opt = gatherUserOptions(opt, varargin{:});
            
            if opt.init;
                K = obj.value.init();
            else
                K = obj.value;
            end
        end
        
        function varargout = assemble(obj, varargin)
            
%                 if obj.assembled
%                     warning('MatrixKernel:assemble', 'The kernel named %s was previously assembled.', obj.name);
%                 end
               
                opt.zero = false;
                opt.tag = obj.options.tag;
%                 opt.component = obj.options.component;
                opt.time = obj.t;
                opt = gatherUserOptions(opt, varargin{:});
                obj.t = opt.time;

                elem = obj.mesh.getElements('tag',opt.tag);
               for i = 1:length(elem);
               
                    if isempty(opt.tag);
                        Ke = obj.evaluateElement(elem(i),obj.t);
                    else
                        Ke = obj.evaluateSide(elem(i), opt.tag, obj.t);
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
        
        function Ke = evaluateElement(obj, elem, t)
              
            if obj.direct            
                Ke = obj.eval(elem, [], t);
                return
            end
            
            if strcmpi(obj.options.type, 'matrix');
                Ke = zeros(elem.n_dof);
            else
                Ke = zeros(elem.n_dof,1);         
            end

            % Loop over the quadrature points
            for i = 1:length(elem.qp);
                Ke = Ke + elem.W(i)*obj.eval(elem, elem.qp{i}, t)*elem.detJ(elem.qp{i});
            end
        end
        
        function Ke = evaluateSide(obj, elem, id, t)
        
            if strcmpi(obj.options.type,'matrix');
                  Ke = obj.evaluateSideMatrix(elem,id,t);
            else
                  Ke = obj.evaluateSideVector(elem,id,t);         
            end
        end
        
        function Ke = evaluateSideMatrix(obj, elem, id, t)
            
            Ke = zeros(elem.n_dof);         

            for s = 1:elem.n_sides; 
                if any(elem.side(s).boundary_id == id);
                    side = elem.buildSide(s);
                    dof = elem.getDof('Side', s, '-local');

                    if elem.local_n_dim == 1;
                        Ke(dof,dof) = Ke(dof,dof) + obj.eval(side,[],t); 
                    else
                        for i = 1:length(side.qp);
                            Ke(dof,dof) = Ke(dof,dof) + side.W(i)*obj.eval(side,side.qp{i},t)*side.detJ(side.qp{i});              
                        end
                    end
                    delete(side);
                end
            end  
        end
        
        function Ke = evaluateSideVector(obj, elem, id, t)

            Ke = zeros(elem.n_dof,1);         
            
            for s = 1:elem.n_sides; 
                if any(elem.side(s).boundary_id == id);
                    side = elem.buildSide(s);
                    dof = elem.getDof('Side', s, '-local');

                    if elem.local_n_dim == 1;
                        Ke(dof) = Ke(dof) + obj.eval(side,[],t); 
                    else
                        for i = 1:length(side.qp);
                            Ke(dof) = Ke(dof) + side.W(i)*obj.eval(side,side.qp{i},t)*side.detJ(side.qp{i});              
                        end
                    end
                    delete(side);
                end
            end  
        end
    end

end