classdef AutoKernel < mFEM.kernels.base.MatrixKernel
    %AUTOMATRIXKERNEL Abstract class for defining finite element matrices

    properties
    end
    
    properties (SetAccess = {?mFEM.registry.MatrixKernelRegistry})
        constReg;
        funcReg;
    end
   
    methods 
        function obj = AutoKernel(mesh, name, eqn, varargin)
            obj = obj@mFEM.kernels.base.MatrixKernel(mesh, name, varargin{:});
            obj.mesh = mesh;
            obj.value = eqn;
            
            obj.options.functions = [];
            obj.options.constants = [];
            [obj.options, unknown] = gather_user_options(obj.options, varargin{:});
            
            if isa(obj.options.constants, 'mFEM.registry.ConstantKernelRegistry');
               obj.constReg = obj.options.constants;
            else
                obj.constReg = mFEM.registry.ConstantKernelRegistry();
            end
            
            if isa(obj.options.functions, 'mFEM.registry.FunctionKernelRegistry');
               obj.funcReg = obj.options.functions;
            else
                obj.funcReg = mFEM.registry.FunctionKernelRegistry();
            end
            
            obj.applyUnknowns(unknown);
            
            obj.parseEquation();
        end
        
        function applyUnknowns(obj, unknown)
           
            for i = 1:2:length(unknown)-1;
               itm = unknown{i};
               val = unknown{i+1};
               
               if isa(itm, 'function_handle');
                   obj.funcReg.add(itm, val, 'constants', obj.constReg);
               elseif ischar(itm);
                   obj.constReg.add(itm, val);
               end
            end
        end               
      
        function parseEquation(obj)
            
            str = obj.value;

            str = regexprep(str,'\<N\>', 'elem.shape(qp)');
            str = regexprep(str,'\<B\>', 'elem.shape_deriv(qp)');      
            if ~strcmp(str, obj.value)
                obj.direct = false; 
                obj.value = str;
            else
                obj.direct = true;
            end

            if obj.direct;
                obj.value = regexprep(obj.value,'\<Ke\>','elem.stiffness()');
            end

            obj.constReg.apply(obj);
            
        end
        
        function value = eval(obj, elem, qp, varargin)
            
            if isempty(varargin); 
                t = [];
            else
                t = varargin{1};
            end
            
            func = obj.value;
            obj.funcReg.apply(func, elem, qp, t);

            value = eval(func);
            
        end
    end

        
        
end