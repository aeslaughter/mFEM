classdef System < mFEM.handle_hide
    %SYSTEM A class for automatic assembly of finite element equations.
    %
    % Syntax
    %   sys = System(mesh)
    %   sys = System(mesh, 'PropertyName', PropertyValue)
    %
    % Description
    %   sys = System(mesh) creates a System object based on the supplied
    %       FEmesh object.
    %   sys = System(mesh, 'PropertyName', PropertyValue) same as above but
    %       allows the options to be changes
    %
    % System Properties Descriptions
    %   Time
    %       true | {false}
    %       A toggle for displaying the matrix and vector assembly time.
    %   
    % See Also FEMESH
    
    properties (SetAccess = private, GetAccess = public)
        mesh = mFEM.FEmesh.empty;   % mesh object
        opt = ...                   % struct of default user options
            struct('time', true);
    end
    
    properties(Access = private)
        reserved = {'N','B','x','y','z','xi','eta','zeta'};
        mat = struct('name', char, 'eqn', char, 'func', char, 'matrix',sparse([]),'boundary_id', uint32([]));
        vec = struct('name', char, 'eqn' ,char, 'func', char, 'vector',[],'boundary_id', uint32([]));
        const = struct('name', char, 'value',[]);
    end
    
    methods (Access = public)
        function obj = System(mesh, varargin)
            %SYSTEM Class constructor.
            
            % Parse the user-defined options
            obj.opt = gather_user_options(obj.opt, varargin{:});
            
            % Store the mesh object
            obj.mesh = mesh;
        end
        
        function add_constant(obj, varargin)
            %ADD_CONSTANT Adds constant(s) variables to the system.
            %
            % Syntax
            %   add_constant('ConstantName', ConstantValue, ...)
            %
            % Description
            %   add_constant('ConstantName', ConstantValue, ...) adds
            %   a single or multiple constants to the System, the
            %   'ConstantName' may be any string identifier and the value
            %   may be a numeric (scalar or matrix) or a string that may be
            %   evaluated using MATLAB's eval function.
            %
            % Examples
            %   sys.add_constant('k',10);
            %   sys.add_cosntant('D','k^2');
            
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
            obj.mat(idx).eqn = eqn;
            obj.mat(idx).func = obj.parse_equation(eqn);
%             obj.mat(idx).matrix = sparse(obj.mesh.n_dof, obj.mesh.n_dof);
            obj.mat(idx).matrix = sparse([]);
            obj.mat(idx).boundary_id = varargin{:};
        end

        function add_vector(obj, name, eqn, varargin)  
            %ADD_VECTOR Create a finite element vector
            
            idx = length(obj.vec) + 1;
            obj.vec(idx).name = name;
            obj.vec(idx).eqn = eqn;
            obj.vec(idx).func = obj.parse_equation(eqn);
            obj.vec(idx).vector = zeros(obj.mesh.n_dof, 1);
            obj.vec(idx).boundary_id = varargin{:};   
        end
        
        function X = get(obj, name)
            %GET Returns the value for the specified name
            
            [type, idx] = obj.locate(name);
            
            switch type;
                case 'matrix'; X = obj.mat(idx).matrix;
                case 'vector'; X = obj.vet(idx).vector;
                case 'constant'; X = obj.const(idx).value;
            end
        end
        
        function X = assemble(obj, name)
            %ASSEMBLE Assembles matrix or vector given by name

            % Locate the matrix or vector
            [type, idx] = obj.locate(name);
            
            % Display a time message
            if obj.opt.time;
                ticID = tmessage('%s %s assembly...', name, type);
            end
            
            % Call the correct assembly routine
            switch lower(type);
               case 'matrix'; X = obj.assemble_matrix(idx); 
               case 'vector'; X = obj.assemble_vector(idx);
               otherwise
                   error('System:Assemble','No assembly routine for %s types', type);
            end
            
            % End message
            if obj.opt.time;
                tmessage(ticID);
            end
        end
    end
    
    methods (Hidden = true, Access = private)
        
        function [type, idx] = locate(obj, name)
            %LOCATE Returns the type and index for the supplied name
            
            % Initialize variables
            idx = [];                       % location of name
            types = {'mat','vec','const'};  % type of entity 
            type_out = {'matrix', 'vector', 'constant'}; % output
            
            % Loop throug the types
            for t = 1:length(types);
                type = type_out{t}; % the current type
                
                % Loop through all array values for the current type
                for i = 1:length(obj.(types{t}));
                    
                    % If name is found return the index
                    if strcmp(obj.(types{t})(i).name, name);
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
            
            % If is string, insert existing constants and evaluate
            if ischar(value);
                value = eval(obj.apply_constants(value));
            end
            
            % Assign the numeric constant
            obj.const(idx).value = value;   % constant value
        end
        
        function fcn = parse_equation(obj, eqn, varargin)
            %PARSE_EQUATION Converts given equation into a useable function
            %
            % Syntax:
            %   parse_equation(eqn)
            %   parse_equation(eqn,'side');
            
            % Get the dimensions of the FE space (use the local, it is what
            % determines the number of xi,eta,... inputs)
            n_dim = obj.mesh.local_n_dim;
            
            % Adjust for special side case
            if nargin == 3 && strcmpi(varargin{1},'side');
                n_dim = n_dim - 1;
            end

            % Build variable strings based on the dimensions of FE space
            if n_dim == 0; % (side of 1D elements)
                var = '';
            elseif n_dim == 1;
                var = 'xi';
            elseif n_dim == 2;
                var = 'xi,eta'; 
            elseif n_dim == 3;
                var = 'xi,eta,zeta';
            else
                error('System:parse_equation', '%d-D finite element space not supported', n_dim);
            end
            
            % Apply constants
            eqn = obj.apply_constants(eqn);
            
            % Insert element shape function and shape function derivatives
            eqn = regexprep(eqn,'N',['elem.shape(',var,')']);
            eqn = regexprep(eqn,'B',['elem.shape_deriv(',var,')']);

            % Create the function string
            if isempty(var);
                fcn = ['@(elem) ', eqn];     
            else
                fcn = ['@(elem,',var,') ', eqn];
            end
        end
        
        function eqn = apply_constants(obj, eqn)
            %APPLY_CONSTANTS applies constants to a string equation
                        
            % Loop through each constant and add if present
            for i = 1:length(obj.const);
                
                % Current constant name and value
                str = obj.const(i).name;    
                val = obj.const(i).value;  

                % Look for the constant as a complete string    
                x1 = regexpi(eqn, str);
                
                % Look for the constant without mathmatical symbols
                s = textscan(eqn, '%s', 'delimiter', '*+-./^'''); s = s{1};
                s = strtrim(s); % remove leading/trailing whitespace
                x2 = find(strcmp(obj.const(i).name, s),1);

                % If the constant is present in both cases above, then
                % insert it into the equation. Using the two cases
                % elimnates problems with having constant names that are
                % within other constant names (e.g., b and Qb)
                if ~isempty(x1) && ~isempty(x2);
                   eqn = regexprep(eqn, str, mat2str(val));
                end
            end
        end
        
        function K = assemble_matrix(obj, idx)
            %ASSEMBLE_VECTOR Assembles a vector for given index

            % Clear existing matrix
            obj.mat(idx).matrix = sparse(obj.mesh.n_dof, obj.mesh.n_dof);
          
            % Case when the vector is only applied to a side
            if ~isempty(obj.mat(idx).boundary_id);
               obj.assemble_matrix_side(idx);
               
            % Case when the vector is applied to entire element    
            else
                obj.assemble_matrix_elem(idx);
            end
            
            % Output the global vector
            K = obj.mat(idx).matrix;
        end   
        
        function assemble_matrix_side(obj, idx)
            %ASSEMBLE_MATRIX_SIDE computes matrix equation on boundary
           
            % Build the function for the side
            fcn = str2func(obj.parse_equation(obj.mat(idx).eqn, 'side'));

            % Extract the boundary ids
            id = obj.mat(idx).boundary_id;
            
            % Create a Matrix object
            matrix = mFEM.Matrix(obj.mesh);
            
            % Loop through each element
            for e = 1:obj.mesh.n_elements;

                % The current element
                elem = obj.mesh.element(e);

                % Intialize the force fector
                Ke = zeros(elem.n_dof, elem.n_dof);
            
                % Loop through all of the boundaries specified
                for i = 1:length(id);

                    % Loop through all sides of this element
                    for s = 1:elem.n_sides; 

                        % Test if this side is on the boundary
                        if any(elem.side(s).boundary_id == id(i));

                            % Create the side element
                            side = elem.build_side(s);

                            % If elem is 1D, then the side is a point that
                            % does not require intergration
                            if elem.local_n_dim == 1;
                                dof = elem.get_dof(s);
                                Ke(dof,dof) = Ke(dof,dof) + fcn(side);

                            % Side elements that are not points    
                            else
                                [qp,W] = side.quad.rules('-cell');

                                % Local dofs for the current side
                                dof = elem.get_dof('Side',s,'-local');

                                % Perform quadrature
                                for j = 1:size(qp,1);
                                    Ke(dof,dof) = Ke(dof,dof) + W(j)*fcn(side,qp{j,:})*side.detJ(qp{j,:});
                                end
                            end

                            % Delete the side element
                            delete(side);
                        end
                    end   
                end 
                     
                % Add the local force vector to the global vector
                dof = elem.get_dof();    
                matrix.add_matrix(Ke, dof, dof);
                %obj.mat(idx).matrix(dof,dof) = obj.mat(idx).matrix(dof,dof) + Ke;
            end
            
            % Build the sparse matrix
            obj.mat(idx).matrix = matrix.init();
        end
        
        function assemble_matrix_elem(obj, idx)
            %ASSEMBLE_MATRIX_ELEM Assemble a matrix for entire domain

            % Create function from the function string
            fcn = str2func(obj.mat(idx).func);

            % Create a Matrix object
            matrix = mFEM.Matrix(obj.mesh.n_dof, obj.mesh.n_dof);

            % Loop through all of the elements
            for e = 1:obj.mesh.n_elements;

                % Extract current element
                elem = obj.mesh.element(e);

                % Initialize the local stiffness matrix
                Ke = zeros(elem.n_dof);

                % Get the quadrature rules for this element
                [qp, W] = elem.quad.rules('cell');

                % Loop through all of the quadrature points and add the
                % result to the local matrix
                for i = 1:size(qp,1);
                    Ke = Ke + W(i)*fcn(elem, qp{i,:})*elem.detJ(qp{i,:});
                end

                % Extract the global dof for this element
                dof = elem.get_dof();   

                % Add the contribution to the global matrix
                matrix.add_matrix(Ke, dof);
                %obj.mat(idx).matrix(dof,dof) = obj.mat(idx).matrix(dof,dof) + Ke;
            end

            % Build the sparse matrix
            obj.mat(idx).matrix = matrix.init();
        end
        
        function f = assemble_vector(obj, idx)
            %ASSEMBLE_VECTOR Assembles a vector for given index
            
            % Clear existing vector
            obj.vec(idx).vector = zeros(obj.mesh.n_dof, 1);
          
            % Case when the vector is only applied to a side
            if ~isempty(obj.vec(idx).boundary_id);
               obj.assemble_vector_side(idx);
               
            % Case when the vector is applied to entire element    
            else
                obj.assemble_vector_elem(idx);
            end
            
            % Output the global vector
            f = obj.vec(idx).vector;
        end  
        
        function assemble_vector_side(obj, idx)
            %ASSEMBLE_VECTOR_SIDE computes vector equation on boundary
            
            % Build the function for the side
            fcn = str2func(obj.parse_equation(obj.vec(idx).eqn, 'side'));
            
            % Extract the boundary ids
            id = obj.vec(idx).boundary_id;
            
            % Loop through each element
            for e = 1:obj.mesh.n_elements;

                % The current element
                elem = obj.mesh.element(e);

                % Loop through all of the boundaries specified
                for i = 1:length(id);

                    % Loop through all sides of this element
                    for s = 1:elem.n_sides; 
                        
                        % Intialize the force fector
                        fe = zeros(elem.n_dof, 1);
                        
                        % Test if this side is on the boundary
                        if any(elem.side(s).boundary_id == id(i));

                            % Create the side element
                            side = elem.build_side(s);

                            % If elem is 1D, then the side is a point that
                            % does not require intergration
                            if elem.local_n_dim == 1;
                                dof = elem.get_dof('Side', s, '-local');
                                fe(dof) = fe(dof) + fcn(side);

                            % Side elements that are not points    
                            else
                                [qp,W] = side.quad.rules('cell');

                                % Local dofs for the current side
                                dof = elem.get_dof('Side',s,'-local');

                                % Perform quadrature
                                for j = 1:size(qp,1);
                                    fe(dof) = fe(dof) + W(j)*fcn(side,qp{j,:})*side.detJ(qp{j,:});
                                end
                            end

                            % Delete the side element
                            delete(side);
                            
                            % Add the local force vector to the global vector
                            dof = elem.get_dof();    
                            obj.vec(idx).vector(dof) = obj.vec(idx).vector(dof) + fe;
                        end
                    end   
                end 
            end
        end
        
        function assemble_vector_elem(obj, idx)
            %ASSEMBLE_VECTOR_ELEM Assebles a vector across entire domain 
            
            % Create the function for element calculations
            fcn = str2func(obj.vec(idx).func);

            % Loop through each element
            for e = 1:obj.mesh.n_elements;

                % The current element
                elem = obj.mesh.element(e);

                % Intialize the force fector
                fe = zeros(elem.n_dof,1);

                % Get the quadrature points in cell form
                [qp, W] = elem.quad.rules('cell');

                % Loop through all of the quadrature points
                for j = 1:size(qp,1);
                    fe = fe + W(j)*fcn(elem,qp{j,:})*elem.detJ(qp{j,:});
                end

                % Add the local force vector to the global vector
                dof = elem.get_dof();    
                obj.vec(idx).vector(dof) = obj.vec(idx).vector(dof) + fe;
            end
        end
    end
end