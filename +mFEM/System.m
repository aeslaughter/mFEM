classdef System < mFEM.handle_hide
    %SYSTEM A class for automatic assembly of finite element equations.
    %
    % Syntax:
    %   sys = System(mesh)
    %
    % Description:
    %
    %

    properties(Access = private)
        mesh = mFEM.FEmesh.empty;
        reserved = {'N','B','x','y','z','xi','eta','zeta'};
        mat = struct('name', char, 'equation', char, 'func',char, 'matrix',sparse([]),'boundary_id', uint32([]));
        vec = struct('name', char, 'equation' ,char, 'func',char, 'vector',[],'boundary_id', uint32([]));
        const = struct('name', char, 'value',[]);
    end
    
    methods 
        function obj = System(mesh)
            %SYSTEM Class constructor.
            obj.mesh = mesh;
        end
        
        function add_constant(obj, varargin)
            %ADD_CONSTANT Adds constant(s) variables to the system.
            %
            % Syntax:
            %   add_constant('ConstantName', ConstantValue, ...)
            
            % Location of last ConstantName
            n = nargin - 2;
            
            % Loop through each name and store in the const property
            for i = 1:2:n;
                obj.add_const(varargin{i}, varargin{i+1});
            end
        end
        
        function add_matrix(obj, name, eqn, varargin)  
            %ADD_MATRIX Create a sparse finite element matrix
            
            % Storate location in matrix array
            idx = length(obj.mat) + 1;
            
            obj.mat(idx).name = name;
            obj.mat(idx).equation = eqn;
            obj.mat(idx).func = obj.parse_equation(eqn);
            obj.mat(idx).matrix = sparse(obj.mesh.n_dof, obj.mesh.n_dof);
            obj.mat(idx).boundary_id = varargin{:};
        end

        function add_vector(obj, name, eqn, varargin)  
            %ADD_VECTOR Create a finite element vector
            
            idx = length(obj.vec) + 1;
            obj.vec(idx).name = name;
            obj.vec(idx).equation = eqn;
            obj.vec(idx).func = obj.parse_equation(eqn);
            obj.vec(idx).vector = zeros(obj.mesh.n_dof, 1);
            obj.vec(idx).boundary_id = varargin{:};   
        end
        
        function X = get(obj, name)
            %GET Returns the value for the specified name
            
            [type, idx] = locate(name);
            
            switch type;
                case 'mat'; X = obj.mat(idx).matrix;
                case 'vec'; X = obj.vet(idx).vector;
                case 'const'; X = obj.const(idx).value;
            end
        end
        
        function X = assemble(obj, name)
            % Assembles matrix and vectors

            % Locate the matrix or vector
            [type,idx] = obj.locate(name);
            
            % Call the correct assembly routine
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
            %LOCATE Returns the type and index for the supplied name
            
            % Initialize variables
            idx = [];                       % location of name
            types = {'mat','vec','const'};  % type of entity 
            
            % Loop throug the types
            for t = 1:length(types);
                type = types{t}; % the current type
                
                % Loop through all array values for the current type
                for i = 1:length(obj.(type));
                    
                    % If name is found return the index
                    if strcmpi(obj.(type)(i).name, name);
                        idx = i;
                        return;
                    end
                end
            end
            
            % Throw and error if the name was not found
            if isempty(idx);
                error('System:locate', 'The entity with name %s was not found.', name);
            end   
        end
        
        function add_const(obj, name, value)
            %ADD_CONST Adds a single constant to the system
            
            % Test that the constant is not reserved
            if any(strcmp(name,obj.reserved));
                error('System:add_const', 'The constant %s is a reserved string, select a different name.', name);
            end
            
            % Add the constant
            idx = length(obj.const) + 1;    % location
            obj.const(idx).name = name;     % constant name
            obj.const(idx).value = value;   % constant value
        end
        
        function fcn = parse_equation(obj, eqn)
            %PARSE_EQUATION Converts given equation into a useable function
            
            % Get the dimensions of the FE space
            n_dim = obj.mesh.n_dim;

            % Build variable strings based on the dimensions of FE space
            if n_dim == 1;
                var = 'xi';
            elseif n_dim == 2;
                var = 'xi,eta'; 
            elseif n_dim == 3;
                var = 'xi,eta,zeta';
            else
                error('System:parse_equation', '%d-D finite element space not supported', n_dim);
            end

            % Insert element shape function and shape function derivatives
            eqn = regexprep(eqn,'N',['elem.shape(',var,')']);
            eqn = regexprep(eqn,'B',['elem.shape_deriv(',var,')']);

            % Loop through each constant and add if present
            for i = 1:length(obj.const);
                
                % Current constant name and value
                str = obj.const(i).name;    
                val = obj.const(i).value;  

                % Look for the constant as a complete string    
                x1 = regexpi(eqn, str);
                
                % Look for the constant without mathmatical symbols
                s = textscan(eqn, '%s', 'delimiter', '*+-./^'''); s = s{1};
                x2 = find(strcmp(obj.const(i).name, s),1);

                % If the constant is present in both cases above, then
                % insert it into the equation. Using the two cases
                % elimnates problems with having constant names that are
                % within other constant names (e.g., b and Qb)
                if ~isempty(x1) && ~isempty(x2);
                   eqn = regexprep(eqn, str, mat2str(val));
                end
            end

            % Create the function string
            fcn = ['@(elem,',var,') ', eqn];
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

        
        
    end
    
 
end