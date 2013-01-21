classdef System < handle
    %SYSTEM A class for automatic assembly of finite element equations.
    % This class allows the specification of the finite element equations
    % for matrices and vectors as strings, the assembly is handled
    % automatically.
    %
    % Example
    %   mesh = FEmesh()
    %   mesh.grid('Line2',0,1,10)
    %   mesh.init()
    %
    %   sys = System(mesh)
    %   sys.add_constant('k',10);
    %   sys.add_matrix('K','B''*k*B');
    %   K = sys.assemble('K');
    %   full(K)
    %
    % See Also FEMESH
    %
    %----------------------------------------------------------------------
    %  mFEM: An Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2012 Andrew E Slaughter
    % 
    %  This program is free software: you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation, either version 3 of the License, or
    %  (at your option) any later version.
    % 
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    % 
    %  You should have received a copy of the GNU General Public License
    %  along with this program.  If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
    
    properties(Access = public);
        time = 0;                   
        
        registry = mFEM.base.Registry;
    end
    
    properties (SetAccess = private, GetAccess = public)
        mesh = mFEM.FEmesh.empty;   % mesh object
        opt = ...                   % struct of default user options
            struct('time', false);
    end
    
    properties(Access = private)
        m_reg;
        c_reg;
        f_reg;
    end

    methods (Access = public)
        function obj = System(mesh, varargin)
            %SYSTEM Class constructor.    
            %
            % Syntax
            %   sys = System(mesh)
            %   sys = System(mesh, 'PropertyName', PropertyValue)
            %
            % Description
            %   sys = System(mesh) creates a System object based on the
            %   supplied FEmesh object.
            %
            %   sys = System(mesh, 'PropertyName', PropertyValue) same as 
            %   above but allows the options to be changes
            %
            % System Properties Descriptions
            %   Time
            %       true | {false}
            %       A toggle for displaying the matrix and vector 
            %       assembly time.
            %   
            % See Also FEMESH
            
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
            %   sys.add_constant('k', 10, 'r', 2);
            %   sys.add_cosntant('D','k^2');
            
            obj.c_reg.add(varargin{:});
        end
        
        function add_function(obj, varargin)
            %ADD_CONSTANT Adds constant(s) variables to the system.
            %
            % Syntax
            %   add_function('FunctionName', FunctionHandle, ...)
            %
            % Description
            %   add_function('FunctionName', FunctionHandle) adds
            %   a functions to the System, the 'FunctionName' may be 
            %   any string identifier and the FunctionHandle may be a
            %   valid MATLAB function handle or character string, see
            %   details below.
            %
            % Description of Function Handle Input
            %   If the FunctionHandle property is a MATLAB function handle
            %   it must be in the following form:
            %
            %       FunctionHandle = @(elem, x, t) ...
            %       
            %   elem = An Element class object.
            %   x = The current position, where x(1) is the x-direction, 
            %   x(2) is the y-direction, and x(3) is the z-direction
            %   t = The current time, which is automatically passed to the
            %   function from the t variable of the System class itself.
            %
            % Description of Text Input
            %   It is also possible to specify a function as a text string
            %   such that the following is valid:
            %       fhandle = str2func(['@(elem,x,t)', <Text Input Here>]);
            %
            % Examples
            %   sys.add_function('k',@(elem,x,t) x(1)^2 + x(2)^2);
            %   sys.add_fucntion('p','x(1)^3');
            
            obj.f_reg.add(varargin{:});
        end

        function add_matrix(obj, varargin)  
            %ADD_MATRIX Create a sparse finite element matrix for assembly
            %
            % Syntax
            %   add_matrix('MatrixName', MatrixEqn, ...)
            %   add_matrix('MatrixName', MatrixEqn, ... 'PropertyName', PropertyValue, ...)
            %
            % Description
            %   add_matrix('MatrixName', MatrixEqn) adds a matrix or
            %   matrices the associated equations for finite element
            %   assembly. The MatrixEqn is a string that is valid MATLAB
            %   (i.e., it works with the eval function) that gives the 
            %   integrands of the finite element equations at the element
            %   level. The variable N is predefined as the element shape 
            %   functions vector and B is pre-defined as grad(N). Any
            %   constant added with the ADD_CONSTANT method may be used.
            %
            %   add_matrix('MatrixName', MatrixEqn, 'PropertyName',
            %       PropertyValue) same as above, but the equation is 
            %   limited according to the supplied options. For example, the
            %   'Boundary' property limits the equation application to the 
            %   boundaries with the same id (see FEMESH.ADD_BOUNDARY)
            %
            % ADD_MATRIX Property Descriptions
            %   Boundary
            %       numeric
            %       Limits the application of the supplied equation to the 
            %       boundaries with the id, see FEMESH.ADD_BOUNDARY
            %
            %   Subdomain
            %       numeric
            %       Limits the application of the supplied equation to the 
            %       elements on the subdomain, see FEMESH.ADD_SUBDOMAIN
            %
            % Examples
            %   sys.add_matrix('M','N''*N');
            %   sys.add_matrix('K','B''*B','Boundary',1);
  
            obj.m_reg.add(varargin{:});
        end

        function add_vector(obj, name, input, varargin)  
            %ADD_VECTOR Create a finite element vector for assembly
            %
            % Syntax
            %   add_vector('VectorName', VectorEqn)
            %   add_vector('VectorName', NumericVector)
            %   add_vector(..., 'PropertyName', PropertyValue)
            %
            % Description
            %   add_vector('VectorName', VectorEqn) adds
            %   a vector and the associated equation for finite element
            %   assembly. The VectorEqn is a string that is valid MATLAB
            %   (i.e., it works with the eval function) that gives the 
            %   integrands of the finite element equations at the element
            %   level. The variable N is predefined as the element shape 
            %   functions vector and B is pre-defined as grad(N). Any
            %   constant added with the ADD_CONSTANT method may be used.
            %
            %   add_vector('VectorName', NumericVector) directly insert a
            %   numeric vector.
            %
            %   add_vector('VectorName', VectorEqn, 'PropertyName',
            %       PropertyValue) same as above, but the equation is 
            %   limited according to the supplied options. For example, the
            %   'Boundary' property limits the equation application to the 
            %   boundaries with the same id (see FEMESH.ADD_BOUNDARY)
            %
            % ADD_VECTOR Property Descriptions
            %   Boundary
            %       numeric
            %       Limits the application of the supplied equation to the 
            %       boundaries with the id, see FEMESH.ADD_BOUNDARY
            %
            %   Subdomain
            %       numeric
            %       Limits the application of the supplied equation to the 
            %       elements on the subdomain, see FEMESH.ADD_SUBDOMAIN
            %
            %   OverWrite
            %       true | {false}
            %       If true the vector is overwritten if the same name is
            %       encounter, otherwise the new vector is added to any
            %       existing.
            %
            %   Dof
            %       numeric vector
            %       Limits the application of the numeric input to specific
            %       dofs.
            %
            % Example
            %   sys.add_vector('f','N''*b'); % b is a constant
            
            % Gather the user options
            options.subdomain = [];
            options.boundary = [];
            options.overwrite = false; 
            options.dof = [];
            options = gather_user_options(options, varargin{:});
            
            % Limit the use of name
            [~,idx] = obj.locate(name, 'mat'); 
            if ~isempty(idx);
                error('ERROR:ADD_VECTOR', 'The name, %s, was previously used for defining a matrix, vectors and matrices may not share names.', name);
            end
            
            % Get the location of vector if overwritting is desired
            if options.overwrite;
                [~,idx] = obj.locate(name,'vec');
            end
            
            % If the index has not been defined, append to any existing
            if isempty(idx);
                idx = length(obj.vec) + 1;
            end
         
            % Append the vector to the structure
            obj.vec(idx).name = name;
            obj.vec(idx).boundary_id = options.boundary; 
            obj.vec(idx).subdomain = options.subdomain;
            
            % Account for different types of vectors
            if isnumeric(input);
               obj.vec(idx).eqn = '';                
                if isempty(options.dof);
                    obj.vec(idx).vector = mFEM.Vector(obj.mesh.n_dof, input); 
                else
                    obj.vec(idx).vector = mFEM.Vector(obj.mesh.n_dof, 0);
                    obj.vec(idx).vector.add_vector(input,options.dof);
                end
            else
               obj.vec(idx).eqn = input; 
               obj.vec(idx).vector = mFEM.Vector.empty;
            end
            
        end
        
        function output = point_value(obj, name, elem, x)
            %POINT_VALUE extract the value for a vector at a point x
            
            [type,idx] = obj.locate(name);
            if ~strcmp(type,'vector');
                error('System:point_value', 'This method only works with vectors');
            end
            
            % Degrees of freedom
            dof = elem.get_dof();

            % Extract local vector
            u = obj.vec(idx).vector.get_local(dof);

            % Shape functions at point x
            N = elem.shape(x{:});

            % Output the value
            output = N*u;
        end
        
        function output = point_gradient(obj, name, elem, x)
            %POINT_GRADIENT extract the gradient of a vector at a point x
            
            [type,idx] = obj.locate(name);
            if ~strcmp(type,'vector');
                error('System:point_value', 'This method only works with vectors');
            end
            
            % Degrees of freedom
            dof = elem.get_dof();

            % Extract local vector
            u = obj.vec(idx).vector.get_local(dof);

            % Shape functions at point x
            B = elem.shape_deriv(x{:});

            % Output the value
            output = B*u;
        end
        
        function X = get(obj, name)
            %GET Returns the value for the specified name.
            %
            % Syntax
            %   X = get(name)
            %
            % Description
            %   X = get(name) returns the value of constant, vector, or 
            %   matrix given the string name. In the case of of matrices
            %   the unassemble Matrix class object is returned.
            
            % Determine the type and location in the structure
            [type, idx] = obj.locate(name);
            
            % Throw and error if the name was not found
            if isempty(idx);
                error('System:locate', 'The entity with name %s was not found.', name);
            end   
            
            % Extract the value based on the type
            switch type;
                case 'matrix';   X = obj.mat(idx).matrix;
                case 'vector';   X = obj.vec(idx).vector;
                case 'constant'; X = obj.const(idx).value;
                case 'function'; X = obj.func(idx).handle;
            end
        end
        
        function [TF, type] = exists(obj, name)
           %EXISTS Returns true if name exists as a type in System
           %
           % Syntax
           %    TF = exists(name)
           %    [TF, type] = exists(name)
           %
           % Description
           %    TF = exists(name) returns a true value if a constant,
           %    matrix, vector, or function exists in the System with the
           %    name given.
           %    [TF, type] = exists(name) same as above but also returns
           %    the type ('mat', 'vec', 'const', or 'func').
            
           % Search for the item
           [type, idx] = obj.locate(name);
           
           % Set the value for TF
           TF = false;
           if ~isempty(idx);
               TF = true;
           end
           
        end
        
        function X = assemble(obj, name)
            %ASSEMBLE Assembles matrix or vector given by name
            %
            % Syntax
            %   X = assemble(name)
            %
            % Description
            %   X = assemble(name) assembles the finite element matrix or
            %   vector. In the case of a matrix the sparse matrix is
            %   returned, for vectors a tradional MATLAB array is given.

            % Locate the matrix or vector
            [type, idx] = obj.locate(name);
            
            % Error if called on something other than a matrix or vector
            if ~any(strcmpi(type,{'matrix','vector'}));
                error('System:Assemble','No assembly routine for %s types', type);
            end
            
            % Display a time message
            if obj.opt.time;
                ticID = tmessage('%s %s assembly...', name, type);
            end
            
            % Call the correct assembly routine
            switch lower(type);
               case 'matrix'; X = obj.assemble_matrix(idx); 
               case 'vector'; X = obj.assemble_vector(idx);
            end
            
            % End message
            if obj.opt.time;
                tmessage(ticID);
            end
        end
    end
    
    methods (Hidden = true, Access = private)
        
        function idx = find_kernel(obj, name)
            
            idx = [];
            for i = 1:length(obj.kernels);
               
                if strcmp(name,obj.kernels{i}.name);
                    idx = i;
                    break;
                end
            end

        end
                
        function [type, idx] = locate(obj, name, varargin)
            %LOCATE Returns the type and index for the supplied name
            % 
            % Syntax
            %   [type,idx] = locate(name)
            %   [type,idx] = locate(name, type)
            %
            % Description
            %   [type,idx] = locate(name) given the string name, the type
            %   ('matrix','vector','constant') and the index within the
            %   storage structure that the name exists.
            %
            %   [type,idx] = locate(name, type) limits the search to the 
            %   specified type, which can be a single type or a cell of
            %   array of types.
            
            % Gather type input, if any
            types = {'mat','vec','const','func'};  % type of entity 
            if nargin == 3;
                if ~iscell(varargin{1});
                    types = varargin(1);
                else
                    types = varargin{1};
                end
            end
            
            % Initialize output index (location of name in the type)
            idx = [];                    
            
            % Loop throug the types
            for t = 1:length(types);
                type = types{t}; % the current type
                
                % Loop through all array values for the current type
                for i = 1:length(obj.(types{t}));
                    
                    % If name is found return the index
                    if strcmp(obj.(types{t})(i).name, name);
                        idx(end+1) = i;
                    end
                end
                
                % It is impossible to have the same name of a different
                % type, so exit if the idx is not empty
                if ~isempty(idx);
                    break;
                end
            end
            
            % Output fullname types
            switch type;
                case 'mat';   type = 'matrix';
                case 'vec';   type = 'vector';
                case 'const'; type = 'constant';
                case 'func';  type = 'function';
            end
        end
        
        function info = test_input_eqn(obj, eqn)
            %TEST_INPUT_EQN
            %

            % Intilize output
            info = struct('direct',{},'function',{},'vector',{},'gradient', {});

            % Search for standard assembly variables
            cstr = {'N','B'};
            for i = 1:length(cstr);                
                [~, S(i)] = obj.insert_string(eqn, cstr{i}, '');
            end
            
            % Search for direct assembly variables
            cstr = {'Ke','fe'};
            for i = 1:length(cstr);
                [~, D(i)] = obj.insert_string(eqn, cstr{i}, '');
            end   
            
            % Test that both are not given
            if any(S) && any(D);
                error('System:test_input_eqn', 'In the equation, %s, both standard assembly (N,B) and direct assembly (Ke,fe) variables were given; it is not possible to mix the assembly types.', eqn);
           
            % Indicate that direct assembly is being used
            elseif any(D);
                info(1).direct = true;
                
            % Indicate that standard assembly is being used
            else
                info(1).direct = false;
            end
            
            % Search for functions
            for i = 1:length(obj.func);
                [~, flag] = obj.insert_string(eqn, obj.func(i).name, '');
                if flag;
                    info(1).function(end+1) = i;
                end 
            end    
            
            % Search for gradient of the vectors
            for i = 1:length(obj.vec);
                expr = ['grad\(\s*',obj.vec(i).name,'\s*\)'];
                ix = regexp(eqn, expr, 'once');
                
                if ~isempty(ix);
                    info(1).gradient(end+1) = i;
                end
            end
            
            % Search for vectors
            for i = 1:length(obj.vec);
                expr = ['(?<!grad\(\s*)', obj.vec(i).name];
                ix = regexp(eqn, expr, 'once');    
                
                if ~isempty(ix);
                    info(1).vector(end+1) = i;
                end 
            end
        end
        
%         function eqn = apply_constants(obj, eqn)
%             %APPLY_CONSTANTS applies constants to a string equation
%             %
%             % Syntax
%             %   eqn = apply_constant(eqn)
%             %
%             % Description
%             %   eqn = apply_constant(eqn) given the eqn string (see 
%             %   ADD_MATRIX and ADD_VECTOR) and applies any constants that
%             %   were given via ADD_CONSTANT.
%                         
%             % Loop through each constant and add if present
%             for i = 1:length(obj.const);
% 
%                 % Current constant name and value
%                 str = obj.const(i).name;    
%                 val = mat2str(obj.const(i).value);  
% 
%                 % Insert the value                
%                 eqn = obj.insert_string(eqn, str, val);
%             end
%         end
        
        function eqn = apply_explicit_vector(obj, eqn, idx, elem, x)
            %APPLY_EXPLICIT_VECTOR Inserts the vector into the equation
            %
            % Syntax
            %   fcn = apply_func(fcn, fidx, elem, x)
            %
            % Description
            %   fcn = apply_func(fcn, fidx, elem, x) applies the function
            %   variables given by the indices in fidx for the position x
            %   to the string fcn supplied.
            
            % Loop through all the vectors
            for i = 1:length(idx);
                
                % The name of the vector to extract from
                name = obj.vec(idx(i)).name;

                % The value extracted from the vector
                value = obj.point_value(name, elem, x);
                
                % Insert the value into the equation
                eqn = obj.insert_string(eqn, name, mat2str(value));             
            end
        end        
        
        function eqn = apply_explicit_vector_grad(obj, eqn, idx, elem, x)
            %APPLY_EXPLICIT_VECTOR Inserts the vector into the equation
            %
            % Syntax
            %   fcn = apply_func(fcn, fidx, elem, x)
            %
            % Description
            %   fcn = apply_func(fcn, fidx, elem, x) applies the function
            %   variables given by the indices in fidx for the position x
            %   to the string fcn supplied.
            
            % Loop through all the vectors
            for i = 1:length(idx);
                
                % The name of the vector to extract from
                name = obj.vec(idx(i)).name;

                % The value extracted from the vector
                value = obj.point_gradient(name, elem, x);
                
                % Insert the value into the equation
                expr = ['grad\(\s*',name,'\s*\)'];
                eqn = regexprep(eqn, expr, mat2str(value));             
            end
        end     
        
        function eqn = apply_func(obj, eqn, idx, elem, x)
            %APPLY_FUNC Inserts the function handle into the equation
            %
            % Syntax
            %   fcn = apply_func(fcn, fidx, elem, x)
            %
            % Description
            %   fcn = apply_func(fcn, fidx, elem, x) applies the function
            %   variables given by the indices in fidx for the position x
            %   to the string fcn supplied.
            
            % Loop through all the functions
            for i = 1:length(idx);
                % Value to insert
                if strcmp(obj.func(i).type,'x');
                    value = feval(obj.func(i).handle, obj, x, obj.time);
                else
                    value = feval(obj.func(i).handle, obj, elem, obj.time);
                end
                
                % Insert the value for the current location
                eqn = obj.insert_string(eqn, ...
                    obj.func(i).name, mat2str(value));             
            end
        end   
        
        function [func, info] = parse_equation(obj, eqn, varargin)
            %PARSE_EQUATION Converts given equation into a useable function
            %
            % Syntax
            %   func = parse_equation(eqn)
            %   func = parse_equation(eqn,'PropertyName',PropertyValue,...)
            %
            % Description
            %   func = parse_equation(eqn) given an equation string, as 
            %   detailed in ADD_MATRIX and ADD_VECTOR, create a working 
            %   equation for use in the assembly routines.
            %
            %   func = parse_equation(eqn,'PropertyName',PropertyValue,...)
            %   same as above but builds the equation according to the
            %   settings described by the properties.
            %
            % PARSE_EQUATION Property Description
            %   side
            %       true | {false}

            % Test the input equation
            info = obj.test_input_eqn(eqn);
            
            % Get the dimensions of the FE space (use the local, it is what
            % determines the number of xi,eta,... inputs)
            n_dim = obj.mesh.local_n_dim;
            
            % Adjust for special side case
            options.side = false;
            options = gather_user_options(options, varargin{:});
            
            % Apply constants
            eqn = obj.apply_constants(eqn);
            
            % Build the traditional shape function equation
            if ~info.direct
            
                %  Adjust the dimension for side case
                if options.side;
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

                % Insert element shape function and shape function derivatives
                eqn = regexprep(eqn,'\<N\>',['elem.shape(',var,')']);
                eqn = regexprep(eqn,'\<B\>',['elem.shape_deriv(',var,')']);
               
            % Build the direct equations    
            else
                eqn = regexprep(eqn,'\<Ke\>','elem.stiffness()');
                eqn = regexprep(eqn,'\<fe\>','elem.force()');
                var = ''; % this causes the proper @ application below
                
            end
            
            % Insert the L variable (L = elem.size())
            eqn = regexprep(eqn,'L','elem.size()');

            % Create the function string
            if isempty(var);
                func = ['@(elem) ', eqn];     
            else
                func = ['@(elem,',var,') ', eqn];
            end
        end

        function K = assemble_matrix(obj, idx)
            %ASSEMBLE_MATRIX Assembles a matrix for given index.
            %
            % Syntax
            %   K = assemble_matrix(idx)
            %
            % Description
            %   K = assemble_matrix(idx) assembles the desired matrix,
            %   where idx is the location in the mat structure.

            % Initialize the output
            K = sparse(obj.mesh.n_dof, obj.mesh.n_dof);
            
            % Loop through the indices
            for i = 1:length(idx);
                
                % Create the Matrix object
                obj.mat(idx(i)).matrix = mFEM.Matrix(obj.mesh);

                % Case when the matrix is only applied to a side
                if ~isempty(obj.mat(idx(i)).boundary_id);
                   obj.assemble_matrix_side(idx(i));

                % Case when the vector is applied to entire element    
                else
                    obj.assemble_matrix_elem(idx(i));
                end

                % Add the matrix to the global matrix
                K = K + obj.mat(idx(i)).matrix.init();
            end
        end   
        
        function assemble_matrix_side(obj, idx)
            %ASSEMBLE_MATRIX_SIDE computes matrix equation on boundary
            %
            % Syntax
            %   K = assemble_matrix_side(idx)
            %
            % Description
            %   K = assemble_matrix_side(idx) assembles the desired matrix 
            %   for a side where idx is the location in the mat structure.
           
            % Parse the equation
            [the_function, info] = obj.parse_equation(obj.mat(idx).eqn,'-side');

            % Extract the boundary ids
            id = obj.mat(idx).boundary_id;
            
            % Get a reference to current Matrix object
            matrix = obj.mat(idx).matrix;
            
            % Loop through each element
            for e = 1:obj.mesh.n_elements;

                % Reset the function string
                str_func = the_function;
                
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
                            
                            % If direct build, perform the build and return
                            % to the next loop iteration
                            if info.direct;
                                dof = elem.get_dof('Side',s,'-local');
                                fcn = str2func(str_func);
                                Ke(dof,dof) = Ke(dof,dof) + fcn(side);
                            
                            % Perform quadrature
                            else
                                % If elem is 1D, then the side is a point that
                                % does not require intergration
                                if elem.local_n_dim == 1;
                                    dof = elem.get_dof(s);
                                    fcn = str2func(str_func);
                                    Ke(dof,dof) = Ke(dof,dof) + fcn(side);

                                % Side elements that are not points    
                                else

                                    % Get side quadrature rule, all sides
                                    if ~exist('qp', 'var')
                                        [qp,W] = side.quad.rules('-cell');
                                    end

                                    % Local dofs for the current side
                                    dof = elem.get_dof('Side',s,'-local');

                                    % Perform quadrature
                                    for j = 1:size(qp,1);

                                        % Apply spatial/temporal functions
                                        if ~isempty(info.function);
                                            x = elem.get_position(qp{i,:});
                                            str_func = obj.apply_func(str_func, info.function, side, x);
                                        end
                
                                        % Apply explicit vectors
                                        if ~isempty(info.vector);
                                            str_func = obj.apply_explicit_vector(str_func, info.vector, side, qp(j,:));
                                        end  
                                        
                                        % Apply explicit vectors gradients
                                        if ~isempty(info.gradient);
                                            str_func = obj.apply_explicit_vector_grad(str_func, info.gradient, side, qp(j,:));
                                        end  
                                        
                                        % Build local matrix
                                        fcn = str2func(str_func);
                                        Ke(dof,dof) = Ke(dof,dof) + W(j)*fcn(side,qp{j,:})*side.detJ(qp{j,:});
                                    end
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
            end
            
            % Store the Matrix class
            obj.mat(idx).matrix = matrix;
        end
        
        function assemble_matrix_elem(obj, idx)
            %ASSEMBLE_MATRIX_ELEM Assemble a matrix for entire domain
            %
            % Syntax
            %   K = assemble_matrix_elem(idx)
            %
            % Description
            %   K = assemble_matrix_elem(idx) assembles the desired matrix 
            %   over the entire domain where idx is the location in the 
            %   mat structure.
            
            % Build the function for the side
            [the_function, info] = obj.parse_equation(obj.mat(idx).eqn);

            % Get a reference to current Matrix object
            matrix = obj.mat(idx).matrix;

            % Get the elements to loop over, if the subdomain it empty it
            % returns all of the elements.
            if isempty(obj.mat(idx).subdomain);
                active_elem = obj.mesh.element;
            else
                active_elem = obj.mesh.get_elements('Subdomain', obj.mat(idx).subdomain);
            end
        
            % Loop through all of the elements
            for e = 1:length(active_elem);
                
                % Reset the function string
                str_func = the_function;

                % Extract current element
                elem = active_elem(e);

                % If direct build, perform the build and return
                % to the next loop iteration
                if info.direct;

                    % Apply spatial/temporal functions
                    if ~isempty(info.function);
                        str_func = obj.apply_func(str_func, info.function, elem, []);
                    end   
                    
                    % Convert to function handle
                    fcn = str2func(str_func);
                    
                    % Compute element matrix
                    Ke = fcn(elem);
                    
                % Otherwise, perform quadrature    
                else
                    % Initialize the local stiffness matrix
                    Ke = zeros(elem.n_dof);

                    % Get the quadrature rules for this and all elements
                    if ~exist('qp','var');
                        [qp,W] = elem.quad.rules('-cell');
                    end

                    % Loop through all of the quadrature points and add the
                    % result to the local matrix
                    for i = 1:size(qp,1);       

                        % Apply spatial/temporal functions
                        if ~isempty(info.function);
                            x = elem.get_position(qp{i,:});
                            str_func = obj.apply_func(str_func, info.function, elem, x);
                        end
                
                        % Apply explicit vectors
                        if ~isempty(info.vector);
                            str_func = obj.apply_explicit_vector(str_func, info.vector, elem, qp(i,:));
                        end  
                        
                        % Apply explicit vectors gradients
                        if ~isempty(info.gradient);
                            str_func = obj.apply_explicit_vector_grad(str_func, info.gradient, elem, qp(i,:));
                        end  
                                        
                        % Case without spatial/temporal functions
                        fcn = str2func(str_func);
                        Ke = Ke + W(i)*fcn(elem, qp{i,:})*elem.detJ(qp{i,:});
                    end
                end
                % Extract the global dof for this element
                dof = elem.get_dof();   

                % Add the contribution to the global matrix
                matrix.add_matrix(Ke, dof);
            end

            % Store the Matrix object
            obj.mat(idx).matrix = matrix;
        end
        
        function f = assemble_vector(obj, idx)
            %ASSEMBLE_VECTOR Assembles a vector for given index
            %
            % Syntax
            %   K = assemble_vector(idx)
            %
            % Description
            %   K = assemble_vector(idx) assembles the desired vector,
            %   where idx is the location in the vec structure.
            
            % Initialize the output
            f = zeros(obj.mesh.n_dof, 1);
            
            % Loop through the indices
            for i = 1:length(idx);
                
                % Skip explicitly defined vectors
                if isempty(obj.vec(idx(i)).eqn);
                    warning('System:asemble_vector', 'There is nothing to assemble for the %s vector', obj.vec(idx(i)).name);
                    f = obj.get(obj.vec(idx(i)).name).init();
                    continue;
                end
                
                % Clear existing vector
                obj.vec(idx(i)).vector = mFEM.Vector(obj.mesh);

                % Case when the vector is only applied to a side
                if ~isempty(obj.vec(idx(i)).boundary_id);
                   obj.assemble_vector_side(idx(i));

                % Case when the vector is applied to entire element    
                else
                    obj.assemble_vector_elem(idx(i));
                end

                % Output the global vector
                f = f + obj.vec(idx(i)).vector.init();
            end
        end  
        
        function assemble_vector_side(obj, idx)
            %ASSEMBLE_VECTOR_SIDE computes vector equation on boundary 
            %
            % Syntax
            %   K = assemble_vector_side(idx)
            %
            % Description
            %   K = assemble_vector_side(idx) assembles the desired vector 
            %   for a side where idx is the location in the vec structure.  
            
            % Parse the equation
            [the_function, info] = obj.parse_equation(obj.vec(idx).eqn,'-side');
            
            % Extract the boundary ids
            id = obj.vec(idx).boundary_id;
            
            % Loop through each element
            for e = 1:obj.mesh.n_elements;

                % Reset the function string
                str_func = the_function;
                
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
                            
                            % Apply explicit vectors
                            if ~isempty(info.vector);
                                dof = elem.get_dof('side',s);
                                str_func = obj.apply_explicit_vector(str_func, info.vector, dof);
                            end
                            
                            % If direct build, perform the build and return
                            % to the next loop iteration
                            if info.direct;
                                dof = elem.get_dof('Side',s, '-local');
                                fcn = str2func(str_func);
                                fe(dof) = fe(dof) + fcn(side);
                                
                            % Perform quadrature    
                            else                            
                                % If elem is 1D, then the side is a point that
                                % does not require intergration
                                if elem.local_n_dim == 1;

                                    % Apply spatial/temporal functions
                                    if ~isempty(info.function);
                                        x = side.nodes;
                                        str_func = obj.apply_func(str_func, info.function, elem, x);
                                    end

                                    % Build local vector
                                    dof = elem.get_dof('Side', s, '-local');
                                    fcn = str2func(str_func);
                                    fe(dof) = fe(dof) + fcn(side);

                                % Side elements that are not points    
                                else
                                    % Get side quadrature rule (non mixed mesh)
                                    if ~exist('qp','var');
                                        [qp,W] = side.quad.rules('-cell');
                                    end

                                    % Local dofs for the current side
                                    dof = elem.get_dof('Side',s,'-local');

                                    % Perform quadrature
                                    for j = 1:size(qp,1);

                                        % Apply spatial/temporal functions
                                        if ~isempty(info.function);
                                            x = elem.get_position(qp{j,:});
                                            str_func = obj.apply_func(str_func, obj.vec(idx).function, side, x);
                                        end
                                        
                                        % Apply explicit vectors
                                        if ~isempty(info.vector);
                                            str_func = obj.apply_explicit_vector(str_func, info.vector, side, qp(j,:));
                                        end
                                        
                                        % Apply explicit vectors gradients
                                        if ~isempty(info.gradient);
                                            str_func = obj.apply_explicit_vector_grad(str_func, info.gradient, side, qp(j,:));
                                        end  
                                        
                                        % Build local vector
                                        fcn = str2func(str_func);
                                        fe(dof) = fe(dof) + W(j)*fcn(side,qp{j,:})*side.detJ(qp{j,:});
                                    end
                                end
                            end
                        
                            % Delete the side element
                            delete(side);
                            
                            % Add the local force vector to the global vector
                            dof = elem.get_dof();    
                            obj.vec(idx).vector.add_vector(fe,dof);
                         end
                    end   
                end 
            end
        end
        
        function assemble_vector_elem(obj, idx)
            %ASSEMBLE_VECTOR_ELEM Assebles a vector across entire domain 
            %
            % Syntax
            %   K = assemble_vector_elem(idx)
            %
            % Description
            %   K = assemble_matrix_elem(idx) assembles the desired matrix 
            %   over the entire domain where idx is the location in the 
            %   vec structure.     
            
            % Build the function for the side
            [the_function, info] = obj.parse_equation(obj.vec(idx).eqn);

            % Get the elements to loop over, if the subdomain is empty it
            % returns all of the elements.
            active_elem = obj.mesh.get_elements('Subdomain', obj.vec(idx).subdomain);

            % Loop through all of the elements
            for e = 1:length(active_elem);

                % Reset the function string
                str_func = the_function;
                
                % Extract current element 
                elem = active_elem(e);
                
                % If direct build, perform the build and return
                % to the next loop iteration
                if info.direct;
                    fcn = str2func(str_func);
                    fe = fcn(side);
                
                % Otherwise perform quadrature    
                else     
                    % Intialize the force fector
                    fe = zeros(elem.n_dof,1);

                    % Get side quadrature rule (non mixed mesh)
                    if ~exist('qp','var');
                        [qp,W] = elem.quad.rules('-cell');
                    end
                    
                    % Loop through all of the quadrature points
                    for j = 1:size(qp,1);

                        % Apply spatial/temporal functions
                        if ~isempty(info.function);
                            x = elem.get_position(qp{j,:});
                            str_func = obj.apply_func(str_func, info.function, elem, x);
                        end
                        
                        % Apply explicit vectors
                        if ~isempty(info.vector);
                            str_func = obj.apply_explicit_vector(str_func, info.vector, elem, qp(j,:));
                        end       
                        
                        % Apply explicit vectors gradients
                        if ~isempty(info.gradient);
                            str_func = obj.apply_explicit_vector_grad(str_func, info.gradient, elem, qp(j,:));
                        end  
                        
                        % Build the local vector
                        fcn = str2func(str_func);
                        fe = fe + W(j)*fcn(elem,qp{j,:})*elem.detJ(qp{j,:});
                    end
                end
                
                % Add the local force vector to the global vector
                dof = elem.get_dof();    
                obj.vec(idx).vector.add_vector(fe,dof);
            end
        end
    end
    
    methods (Static, Hidden = true)
        function [out, found] = insert_string(eqn, str, val)
            %INSERT_STRING Inserts the variables into eqn string
            %
            % Syntax
            %   [eqn, found] = insert_string(eqn, str, val)
            %
            % Description
            %   [eqn, found] = insert_string(eqn, str, val) searchs the eqn
            %   string for str and replaces it with val, if it is found. If
            %   the str is located, then the found output variable is set
            %   to true.

            % Look for the constant as a complete string    
            out = regexprep(eqn, ['\<',str,'\>'], val);

            % Test if output is different
            found = ~strcmp(eqn, out);
        end
    end
end