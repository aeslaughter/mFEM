classdef AutoKernel < mFEM.kernels.base.MatrixKernel
    %AUTOMATRIXKERNEL Abstract class for defining finite element matrices
    
    properties (GetAccess = public, SetAccess = protected) %(SetAccess = {?mFEM.registry.MatrixKernelRegistry})
        const;
        func;
        input;
        vec;
        eqn;
    end
   
    methods 
        function obj = AutoKernel(mesh, name, eqn_input, varargin)
            obj = obj@mFEM.kernels.base.MatrixKernel(mesh,name,varargin{:});
            obj.mesh = mesh;
            obj.input = eqn_input;
            
            obj.options.funcregistry = [];
            obj.options.constantregistry = [];
            [obj.options, unknown] = gatherUserOptions(obj.options, varargin{:});
            
            if isa(obj.options.constantregistry, 'mFEM.registry.ConstantRegistry');
               obj.const = obj.options.constantregistry;
            else
                obj.const = mFEM.registry.ConstantRegistry();
            end
            
            if isa(obj.options.funcregistry, 'mFEM.registry.FuncRegistry');
               obj.func = obj.options.funcregistry;
            else
                obj.func = mFEM.registry.FuncRegistry();
            end
            
            obj.applyUnknowns(unknown);
            
            obj.parseEquation();
        end
        
        function applyUnknowns(obj, unknown)
           
            for i = 1:2:length(unknown)-1;
               itm = unknown{i};
               val = unknown{i+1};
               
               if isa(itm, 'function_handle');
                   obj.func.add(itm, val, 'ConstantRegistry', obj.const);
               elseif ischar(itm);
                   obj.const.add(itm, val);
               end
            end
        end               
      
        function parseEquation(obj)
            
            str = obj.input;
            str = regexprep(str,'\<N\>', 'elem.shape(qp)');
            str = regexprep(str,'\<B\>', 'elem.shapeDeriv(qp)');    
            
            if ~strcmp(str, obj.input)
                obj.direct = false; 
            else
                obj.direct = true;
            end

            if obj.direct;
                str = regexprep(str,'\<Ke\>','elem.stiffness()');
            end

            str = regexprep(str,'\<L\>', 'elem.size()');

            obj.eqn = obj.const.apply(str);
            
        end
        
        function value = eval(obj, elem, qp, t)
            str = obj.func.apply(obj.eqn, elem, qp, t);            
            value = eval(str);
        end
    end   
end