classdef MatrixKernel < mFEM.kernels.base.Kernel                     
    %MATRIXKERNEL Abstract class for defining finite element matrices

    properties
       t = 0;
       options = struct(...
           'tag',[],'component',[],'type','matrix','parallel',false);
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

            [obj.options, ~] = gatherUserOptions(obj.options, varargin{:});
            obj.mesh = mesh;
            
            if any(strcmpi(obj.options.type,{'matrix','mat','m'}));
                obj.options.type = 'matrix';
            elseif any(strcmpi(obj.options.type,{'vector','vec','v'}));
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
        
        function varargout = assemble(obj,varargin)
            
%                 if obj.assembled
%                     warning('MatrixKernel:assemble', 'The kernel named %s was previously assembled.', obj.name);
%                 end

            opt.tag = obj.options.tag;
            opt.parallel = obj.options.parallel;
            opt.zero = false;
            opt.time = obj.t;               
            opt = gatherUserOptions(opt,varargin{:});
            obj.t = opt.time;
            type = obj.options.type;

            int_type = 'element';
            if ~isempty(opt.tag);
                tag = obj.mesh.getTag(opt.tag);
                if isempty(tag);
                    error('MatrixKernel:assemble:IvalidTag','The tag %s was not found.',opt.tag);
                end
                if strcmpi(tag.type,'boundary');
                    int_type = 'side';
                end
            end
                  
            elem = obj.mesh.getElements('tag',opt.tag,'-parallel');
            spmd
                n_dof = sum([elem.n_dof]);
                if strcmpi(type,'vector');
                    K = mFEM.Vector(n_dof);
                elseif strcmpi(type,'matrix');
                    K = mFEM.Matrix(n_dof);
                end

                if strcmpi(int_type,'element');
                    for i = 1:length(elem);
                        Ke = obj.evaluateElement(elem(i),obj.t);
                        dof = elem(i).getDof();
                        K.add(Ke,dof); 
                    end
                    
                elseif strcmpi(int_type,'side');
                     for i = 1:length(elem);
                        Ke = obj.evaluateSide(elem(i), opt.tag, obj.t);
                        dof = elem(i).getDof();
                        K.add(Ke,dof); 
                     end
                end
            end

            obj.value = mFEM.pMatrix(K);

            if nargout == 1;
                varargout{1} = obj.value.assemble('parallel',opt.parallel); 
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
            [~,sid] = elem.hasTag(id);
            for i = 1:length(sid); 
                s = sid(i);
                side = elem.buildSide(s);
                dof = elem.getDof('Side', s, '-local');

                if elem.n_dim == 1;
                    Ke(dof,dof) = Ke(dof,dof) + obj.eval(side,[],t); 
                else
                    for j = 1:length(side.qp);
                        Ke(dof,dof) = Ke(dof,dof) + side.W(j)*obj.eval(side,side.qp{j},t)*side.detJ(side.qp{j});              
                    end
                end
                delete(side);
            end  
        end
        
        function Ke = evaluateSideVector(obj, elem, id, t)

            Ke = zeros(elem.n_dof,1);         
            
            [~,sid] = elem.hasTag(id);
            for i = 1:length(sid); 
                s = sid(i);
                side = elem.buildSide(s);
                dof = elem.getDof('Side', s, '-local');

                if elem.n_dim == 1;
                    Ke(dof) = Ke(dof) + obj.eval(side,[],t); 
                else
                    for j = 1:length(side.qp);
                        Ke(dof) = Ke(dof) + side.W(j)*obj.eval(side,side.qp{j},t)*side.detJ(side.qp{j});              
                    end
                end
                delete(side);
            end 
        
        end
    end

end