classdef System < handle
   
    properties (Access = public)
       time = []; 
    end
    
    properties
        mat = struct([]);;%struct('name', char, 'equation', char, 'func',char, 'matrix',sparse([]));
        vec;% = struct('name', char, 'equation' ,char, 'func',char, 'vector',[]);
        const;% = struct('name', char, 'value',[]);
    end
    
    properties(Access = private)
       mesh; 
    end
    
    methods 
        function obj = System(mesh)
            obj.mesh = mesh;
        end
        
        function add_constant(obj, name, value)
            idx = length(obj.mat) + 1;
            obj.const(idx).name = name;
            obj.const(idx).value = value;
        end
        
        function add_matrix(obj, name, eqn)     
            idx = length(obj.mat) + 1;
            obj.mat(idx).name = name;
            obj.mat(idx).equation = eqn;
            obj.mat(idx).func = obj.parse_equation(eqn);
            obj.mat(idx).matrix = sparse(obj.mesh.n_dof,obj.mesh.n_dof);   
        end
        
        function M = get_matrix(obj, name)
            for i = 1:length(obj.mat);
                if strcmpi(name, obj.mat(i).name);
                    M = obj.mat(i).matrix;
                    return;
                end
            end
        end
        
        function assemble(obj)
           obj.assemble_matrix; 
        end
 
    end
    
    methods (Access = private)
        function func = parse_equation(obj, eqn)
           n_dim = obj.mesh.n_dim;
           
           if n_dim == 1;
               
           elseif n_dim == 2;
               var = 'xi,eta';
                   
           elseif n_dim == 3;

           else

           end

           eqn = regexprep(eqn,'N',['N(',var,')']);
           eqn  = regexprep(eqn,'B',['B(',var,')']);
         
           % constants
         % mat2str(5*eye(2))
         
         func = ['@(',var,') ', eqn];
           
        end
        
        
        
        function assemble_matrix(obj)
           
            qrule = Gauss(2);
            [qp,W] = qrule.rules();
            
            for m = 1:length(obj.mat);
                for e = 1:obj.mesh.n_elements;

                    elem = obj.mesh.element(e);

                    N = @(xi,eta) elem.shape(xi,eta);
                    B = @(xi,eta) elem.shape_deriv(xi,eta);
                    M = str2func(obj.mat(m).func)

                    K = zeros(elem.n_dof);
                    for i = 1:length(qp);
                        for j = 1:length(qp);
                            K = K + W(j)*W(i)*B(qp(i),qp(j))'*B(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
                        end
                    end
                    
                    dof = elem.global_dof;    
                    obj.mat(m).matrix(dof,dof) = K;
                end
               
            end
 
        end
        
        
    end
    
 
end