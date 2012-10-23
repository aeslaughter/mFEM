classdef System < handle
   
    properties (Access = public)
       time = []; 
    end
    
    properties(Access = private)
       mesh; 
        mat = struct('name', char, 'equation', char, 'func',char, 'matrix',sparse([]),'boundary_id', uint32([]));
        vec = struct('name', char, 'equation' ,char, 'func',char, 'vector',[],'boundary_id', uint32([]));
        const = struct('name', char, 'value',[]);
    end
    
    methods 
        function obj = System(mesh)
            obj.mesh = mesh;
        end
        
        function add_constant(obj, varargin)
            
            n = nargin - 1;
            for i = 1:2:n-1;
               obj.add_const(varargin{i},varargin{i+1}); 
            end
        end
        
        function add_matrix(obj, name, eqn, varargin)  
            
            idx = length(obj.mat) + 1;
            
            obj.mat(idx).name = name;
            obj.mat(idx).equation = eqn;
            obj.mat(idx).func = obj.parse_equation(eqn);
            obj.mat(idx).matrix = sparse(obj.mesh.n_dof, obj.mesh.n_dof);
            obj.mat(idx).boundary_id = varargin{:};
        end

        function add_vector(obj, name, eqn, varargin)  
            
            idx = length(obj.vec) + 1;
            
            obj.vec(idx).name = name;
            obj.vec(idx).equation = eqn;
            obj.vec(idx).func = obj.parse_equation(eqn);
            obj.vec(idx).vector = zeros(obj.mesh.n_dof, 1);
            obj.vec(idx).boundary_id = varargin{:};
            
        end
        
%         function M = get_matrix(obj, name)
%             for i = 1:length(obj.mat);
%                 if strcmpi(name, obj.mat(i).name);
%                     M = obj.mat(i).matrix;
%                     return;
%                 end
%             end
%         end
        
        function X = assemble(obj, name)
            

           [type,idx] = obj.locate(name);

           
           switch lower(type);
               case 'mat'; X = obj.assemble_matrix(idx); 
               case 'vec'; X = obj.assemble_vector(idx);
               otherwise
                   error('System:Assemble','No assembly routine for %s types', type);
           end
        end
 
    end
    
    methods (Access = private)
        
        function [type, idx] = locate(obj, name)
            
            
            idx = [];
            types = {'mat','vec','const'};
            for t = 1:length(types);
                type = types{t};
                
                for i = 1:length(obj.(type));
                    if strcmpi(obj.(type)(i).name, name);
                        idx = i;
                        return;
                    end
                end
            end
                 
                 
        end
        
        
        function fcn = parse_equation(obj, eqn)

           n_dim = obj.mesh.n_dim;
           
           if n_dim == 1;
               var = 'xi';
           elseif n_dim == 2;
               var = 'xi,eta'; 
           elseif n_dim == 3;
                var = 'xi,eta,zeta';
           else
                error('System:parse_equation', '%d-D finite element space not supported', n_dim);
           end
           
           eqn = regexprep(eqn,'N',['elem.shape(',var,')']);
           eqn  = regexprep(eqn,'B',['elem.shape_deriv(',var,')']);
           
            for i = 1:length(obj.const);
               str = obj.const(i).name;
               val = obj.const(i).value;
                
                
               x1 = regexpi(eqn, str);
               s = textscan(eqn, '%s', 'delimiter', '*+-./^'''); s = s{1};
               x2 = find(strcmp(obj.const(i).name, s));
               
               if ~isempty(x1) && ~isempty(x2);
                   eqn = regexprep(eqn, str, mat2str(val));
               end

                 
            end
           
           fcn = ['@(elem,',var,') ', eqn]
        end
        
        
        
        
        function K = assemble_matrix(obj, idx)

                fcn = str2func(obj.mat(idx).func);
            
                for e = 1:obj.mesh.n_elements;

                    elem = obj.mesh.element(e);
                    
                    [qp, W] = elem.quad.rules();
                    
                    
                    Ke = zeros(elem.n_dof);
                    for i = 1:length(qp);
                        Ke = Ke + W(i)*fcn(elem,qp(i))*elem.detJ(qp(i));
                    end
                    
                    dof = elem.get_dof();    
                    obj.mat(idx).matrix(dof,dof) = obj.mat(idx).matrix(dof,dof) + Ke;
                end
                
                K = obj.mat(idx).matrix;
                
        end
        
        function f = assemble_vector(obj, idx)

                V = obj.vec(idx);
                fcn = str2func(V.func);
            
                for e = 1:obj.mesh.n_elements;

                    elem = obj.mesh.element(e);
                    
                    [qp, W] = elem.quad.rules();
                                            
                    fe = zeros(elem.n_dof,1);

                    if ~isempty(V.boundary_id);
                        
                        for i = 1:length(V.boundary_id);
                            id = V.boundary_id(i);
                        
                            for s = 1:elem.n_sides; 
                                if elem.side(s).boundary_id == id;
                                    xi = elem.lims(elem.get_dof(s));    % xi value to evaluate at
                                    fe = fe + fcn(elem,xi);              
                                end
                            end   
                        end
                        
                    else
                        for i = 1:length(qp);
                            fe = fe + W(i)*fcn(elem,qp(i))*elem.detJ(qp(i));
                        end
                    end
                    
                    dof = elem.get_dof();    
                    obj.vec(idx).vector(dof) = obj.vec(idx).vector(dof) + fe;
                end
                
                f = obj.vec(idx).vector;
                
        end
        
        function add_const(obj, name, value)
            idx = length(obj.const) + 1;
            obj.const(idx).name = name;
            obj.const(idx).value = value;
        end
        
        
    end
    
 
end